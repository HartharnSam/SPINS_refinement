#!/bin/bash
#SBATCH --time=47:59:00
#SBATCH --array=1%1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
module load intel Autoconf/2.69-GCCcore-8.2.0 gettext/0.19.8.1-GCCcore-8.2.0 Automake/1.16.1-GCCcore-8.2.0 libtool/2.4.6-GCCcore-8.2.0

#Set Directory
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" case_list)
cd $DIR

#run the script
srun ./ice_covered_water.x;
srun ./derivatives.x --config=spins.conf_deriv;
#Move stuff to the folder I want it in
mkdir transfer
mv rho.* ./transfer/;
mv u.* ./transfer/;
mv w.* ./transfer/;
mv diss.* ./transfer/;
mv vorty.* ./transfer/;
mv diagnostics.txt ./transfer/;
mv spins.conf ./transfer/;
mv stresses_bottom* ./transfer/;
mv stresses_top* ./transfer/;
mv *grid* ./transfer/;
mv plot* ./transfer/;
mv tracer* ./transfer/;

cd /nobackup/b5006861/SPINS/SPINS_refinement;
