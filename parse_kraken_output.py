#! /Users/jenniferlspillane/opt/anaconda3/bin/python

import argparse


def calculate_percents(kmer_sets):
    #setting up the lists of both kinds of kmers with 0s so that it's never empty for the math later
    unclassified_kmers = [0.0]
    classified_kmers = [0.0]
    percentage = 0.0
    switch = True
    for item in kmer_sets:
        #split on the colons to separate the id from the number
        num_kmers = item.split(":")
        #if the id is 0 (unclassified), we add that number of kmers to the list, if it's an "A", we ignore it, every other kmer number gets added to the other list.
        if num_kmers[0].strip() == "0":
            unclassified_kmers.append(float(num_kmers[1]))
        elif num_kmers[0].strip() == "A":
            continue
        else:
            classified_kmers.append(float(num_kmers[1]))
            
    total_kmers = ((sum(classified_kmers)) + (sum(unclassified_kmers)))
    #print(total_kmers)
    if total_kmers == 0.0:
        switch = False
    else:       
        #calculate the percentage of classified kmers in the total kmers from the read for the first read
        percentage = (sum(classified_kmers)) / total_kmers
    #print(percentage, switch)
    return [percentage, switch]



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
                
                        #splits into the kmer assignments for each read in the pair
                        assignment = fields[4].split("|:|")
                        #need to split again on whitespace to get the individual id and number pairs
                        kmer_sets1 = (assignment[0].strip()).split()
                        kmer_sets2 = (assignment[1].strip()).split()
                        
                        #call the calculation function for each read in the pair and store the result
                        percent_classified1 = calculate_percents(kmer_sets1)
                        percent_classified2 = calculate_percents(kmer_sets2)

                        if percent_classified1[1] and percent_classified2[1]:
                            #make a tuple to keep everything together and append to the list of all the lines we'll write to the output file
                            lines_for_output.append((read_class, read_name, percent_classified1[0], percent_classified2[0]))
                        else:
                            continue

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
                        
                        #splits into the kmer assignments for each read in the pair
                        kmer_sets3 = fields[4].split()
                        
                        #call the calculation function for the read and store the result
                        percent_classified = calculate_percents(kmer_sets3)
                        
                        if percent_classified[1]:
                            #make a tuple to keep everything together and append to the list of all the lines we'll write to the output file
                            lines_for_output.append((read_class, read_name, percent_classified[0]))
                        else:
                            continue

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
        
        
        

























