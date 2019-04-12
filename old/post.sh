#!/bin/bash
# bash script for submitting a Matlab job to the sharcnet Graham queue

#SBATCH --mem-per-cpu=8G      # memory per processor (default in Mb)
#SBATCH --time=00-01:00         # time (DD-HH:MM)
#SBATCH --job-name="post_iswlong"     # job name
#SBATCH --input=post_iswlong.m   # Matlab script
##SBATCH —dependency=afterok:<jobid>  # Wait for job to complete

#SBATCH --ntasks=1              # number of processors
##SBATCH --nodes=1               # number of nodes
##SBATCH --ntasks-per-node=32    # processors per node
#SBATCH --output=post-%j.log                 # log file
#SBATCH --error=post-%j.err                  # error file
#SBATCH --mail-user=c2xu@uwaterloo.ca     # who to email
#SBATCH --mail-type=FAIL                    # when to email
#SBATCH --account=ctb-mmstastn              # UW Fluids designated resource allocation

module load matlab/2017a
matlab -nodisplay -nosplash -singleCompThread
