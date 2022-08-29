#! /Users/jenniferlspillane/opt/anaconda3/bin/python

import argparse


def parse_kraken():
    lines_for_output = []
    try:
        with open(args.kraken_output, "r") as kraken_file:
            with open(args.output, "w") as output_file:
                #an option that stores "true", so that the function can handle paired or unpaired reads
                if args.paired == True:
            
                    for line in kraken_file:
                        #splitting the lines on tabs first
                        fields = line.split('\t')
                        #this stores the letter at the beginning of the line, either "C" or "U"
                        read_class = fields[0].strip()
                        #the read name is harvested from the header by kraken2 - should be useful for pulling reads later
                        read_name = fields[1].strip()
                
                
                        #setting up some empty strings to hold percentages of classified kmers (because I didn't know how to set up empty floats)
                        percentage1 = ""
                        percentage2 = ""
                        #splits into the kmer assignments for each read in the pair
                        assignment = fields[4].split("|:|")
                        #need to split again on whitespace to get the individual id and number pairs
                        kmer_sets1 = (assignment[0].strip()).split()
                        kmer_sets2 = (assignment[1].strip()).split()
                        #setting up the lists of both kinds of kmers (for both reads in the pair) with 0s so that it's never empty for the math later
                        unclassified_kmers1 = [0.0]
                        classified_kmers1 = [0.0]
                        unclassified_kmers2 = [0.0]
                        classified_kmers2 = [0.0]
                    
                        for item in kmer_sets1:
                            #split AGAIN on the colons to separate the id from the number
                            num_kmers = item.split(":")
                            #if the id is 0 (unclassified), we add that number of kmers to the list, if it's an "A", we ignore it, every other kmer number gets added to the other list.
                            if num_kmers[0].strip() == "0":
                                unclassified_kmers1.append(float(num_kmers[1]))
                            elif num_kmers[0].strip() == "A":
                                continue
                            else:
                                classified_kmers1.append(float(num_kmers[1]))
                    
                        #calculate the percentage of classified kmers in the total kmers from the read for the first read
                        percentage1 = (sum(classified_kmers1)) / ((sum(classified_kmers1)) + (sum(unclassified_kmers1)))

                        #all the same stuff, but for the second read in the pair
                        for item in kmer_sets2:
                            num_kmers = item.split(":")
                            if num_kmers[0].strip() == "0":
                                unclassified_kmers2.append(float(num_kmers[1]))
                            elif num_kmers[0].strip() == "A":
                                continue
                            else:
                                classified_kmers2.append(float(num_kmers[1]))
                            
                        percentage2 = (sum(classified_kmers2)) / ((sum(classified_kmers2)) + (sum(unclassified_kmers2)))

                        #make a tuple to keep everything together and append to the list of all the lines we'll write to the output file
                        lines_for_output.append((read_class, read_name, percentage1, percentage2))

                    #write the output file so that the format works for ggplot (each read in the pair on a different line)
                    for read_line in lines_for_output:
                        output_file.write("{0}\t{1}\tread1\t{2}\n{0}\t{1}\tread2\t{3}\n".format(read_line[0], read_line[1], read_line[2], read_line[3]))

                #now do the same thing, but without the paired reads, the format of the file is different
                else:
                    for line in kraken_file:
                        #splitting the lines on tabs first
                        fields = line.split('\t')
                        #this stores the letter at the beginning of the line, either "C" or "U"
                        read_class = fields[0].strip()
                        #the read name is harvested from the header by kraken2 - should be useful for pulling reads later
                        read_name = fields[1].strip()
                        percentage = ""
                        kmer_sets = fields[4].split()
                        unclassified_kmers = [0.0]
                        classified_kmers = [0.0]
                        for item in kmer_sets:
                            num_kmers = item.split(":")
                            if num_kmers[0].strip() == "0":
                                unclassified_kmers.append(float(num_kmers[1]))
                            elif num_kmers[0].strip() == "A":
                                continue
                            else:
                                classified_kmers.append(float(num_kmers[1]))
                            
                        percentage = (sum(classified_kmers)) / ((sum(classified_kmers)) + (sum(unclassified_kmers)))

                        lines_for_output.append((read_class, read_name, percentage))

                    #output file can be more basic as there aren't as many things to keep straight
                    for read_line in lines_for_output:
                        output_file.write("{0}\t{1}\t{2}\n".format(read_line[0], read_line[1], read_line[2]))



    except IOError:
        print("Issue reading or writing file")


parser = argparse.ArgumentParser(description = "Arguments for parsing the kraken2 stdout output")
parser.add_argument("-k", "--kraken_output", help = "path to the file containing the kraken output for individual reads/sets")
parser.add_argument("-p", "--paired", action = "store_true", help = "a flag to signal that the kraken output comes from paired-end reads")
parser.add_argument("-o", "--output", help = "path to a csv file containing the ratios of classified and unclassified kmers for each read/set")
args = parser.parse_args()

parse_kraken()
                    
#./parse_kraken_output.py -k top_test_conus_kraken.out -p -o test_percentages1.csv
#./parse_kraken_output.py -k top_test_pleuro_kraken.out -o test_percentages2.csv
        
        
        

























