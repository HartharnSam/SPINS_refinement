#!/bin/bash
# Sample bash script for submitting a Matlab job to the sharcnet Graham queue

#SBATCH --mem-per-cpu=2G                         # memory per processor
#SBATCH --time=00-01:00                          # time (DD-HH:MM)
#SBATCH --job-name="matlab"                      # job name
#SBATCH --input=driver_mode1_mode2.m             # Matlab script

#SBATCH --ntasks=1                               # number of processors
#SBATCH --output=logs/mat-%j.log                 # log file
#SBATCH --error=logs/mat-%j.err                  # error file
# #SBATCH --mail-user=username@uwaterloo.ca      # who to email
# #SBATCH --mail-type=FAIL                       # when to email
#SBATCH --account=ctb-mmstastn                   # UW Fluids designated resource allocation

module load matlab/2017a 
matlab -nodisplay -nosplash -singleCompThread


## Add your casenames here. Must match the casenames from the .txt cases file
casenames='base mode1 mode2 gravity lamp'

home=`pwd` #Remember current directory
for casename in $casenames
do
    cp wave_reader.x ../$casename
    cp submit.sh ../$casename

    cd ../$casename
    sbatch submit.sh
    cd $home
done

