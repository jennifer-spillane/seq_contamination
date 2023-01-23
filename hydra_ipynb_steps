### To Get on a Jupyter Notebook 

I mostly followed the instructions for getting on one on Hydra (https://confluence.si.edu/display/HPC/Running+Jupyter+on+Hydra), but I'm just going to record it here too so that I have a streamlined version for this new machine I'm working on. I also put these instructions into a notebook on this machine, but that's not a very good solution for actually using it, so, I'll duplicate here too.  

I need to make sure I have an environment activated that has Jupyter Lab, obviously, or this whole thing won't work.
`mamba activate ml_env`  

Then I run this command to launch the thing, although if I already have a jupyter notebook running locally, I should change the port to something else.    
`jupyter lab --no-browser --ip=`hostname` --port=8888`  *Note that this will not run as is, because markdown takes the necessary puctuation marks away. Better to double click and use the commands from inside the cell.*  

Next, in a terminal tab on my local computer, I run this command and enter my password when prompted.  
`ssh -N -L 8888:boxx:8888 spillanej@hydra-dl04.si.edu`  

And finally, I go to this site: http://localhost:8888. It will prompt me for a token, and I enter the token that the notebook gave me in the boxx window when I ran the jupyterlab command. It should be in a line like this `http://boxx:8888/lab?token=235a1e53f67b8a3b46eb74a9a82eb92b3b8f910d2ab316df`, so just the last part is what I need. I could set up a password too, but honestly it's not that big of a deal to enter this, so that's how I'm going to do it for now.
