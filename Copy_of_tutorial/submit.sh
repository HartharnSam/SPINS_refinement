#!/bin/bash
# bash script that matlab.sh uses to submit a job to SLURM for each cases

#SBATCH --mem-per-cpu=2G                    # memory per processor (default in Mb)
#SBATCH --time=00-01:00                     # time (DD-HH:MM)
#SBATCH --job-name="mode1-2"                # job name

#SBATCH --ntasks=4                          # number of processors
#SBATCH --output=sim-%j.log                 # log file
#SBATCH --error=sim-%j.err                  # error file
# #SBATCH --mail-user=user@uwaterloo.ca     # who to email
# #SBATCH --mail-type=FAIL                  # when to email
#SBATCH --account=ctb-mmstastn              # UW Fluids designated resource allocation

srun ./wave_reader.x

