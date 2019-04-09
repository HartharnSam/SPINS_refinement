#!/bin/bash
# bash script for submitting a Matlab job to the sharcnet Graham queue

#SBATCH --mem-per-cpu=2G      # memory per processor (default in Mb)
#SBATCH --time=00-01:00         # time (DD-HH:MM)
#SBATCH --job-name="matlab"     # job name
#SBATCH --input=driver_tilttank.m  # Matlab script
##SBATCH —dependency=afterok:<jobid>  # Wait for job to complete

#SBATCH --ntasks=1              # number of processors
##SBATCH --nodes=1               # number of nodes
##SBATCH --ntasks-per-node=32    # processors per node
#SBATCH --output=../dummy/mat-%j.log                 # log file
#SBATCH --error=../dummy/mat-%j.err                  # error file
#SBATCH --mail-user=a2grace@uwaterloo.ca     # who to email
#SBATCH --mail-type=FAIL                    # when to email
#SBATCH --account=ctb-mmstastn              # UW Fluids designated resource allocation

module load matlab/2017a 
matlab -nodisplay -nosplash -singleCompThread
 
casenames='ssamp smamp lmamp llamp'

for casename in $casenames
do
    cp wave_reader.x ../$casename
    cp submit.sh ../$casename
#   cp post_iswlong.m ../$casename
#   cp post.sh ../$casename

    cd ../$casename
    sbatch submit.sh
    cd ../matlab2spins
done

