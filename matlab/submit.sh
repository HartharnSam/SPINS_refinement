#!/bin/bash
# bash script for submitting a Matlab job to the sharcnet Graham queue

#SBATCH --mem-per-cpu=2G    # memory per processor (default in Mb)
#SBATCH --time=7-00:00         # time (DD-HH:MM)
#SBATCH --job-name="ttank"     # job name
##SBATCH --input=post_iswlong.m   # Matlab script
##SBATCH --dependency=afterok:<jobid>  # Wait for job to complete

#SBATCH --ntasks=256            # number of processors
##SBATCH --nodes=1               # number of nodes
##SBATCH --ntasks-per-node=32    # processors per node
#SBATCH --output=sim-%j.log                 # log file
#SBATCH --error=sim-%j.err                  # error file
#SBATCH --mail-user=a2grace@uwaterloo.ca     # who to email
#SBATCH --mail-type=FAIL                    # when to email
#SBATCH --account=ctb-mmstastn              # UW Fluids designated resource allocation

srun ./wave_reader.x

#sbatch post.sh

