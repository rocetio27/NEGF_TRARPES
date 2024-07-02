#!/bin/bash
#SBATCH -J tt         # job name
#SBATCH -o myMPI.o%j       # output and error file name (%j expands to jobID)
#SBATCH -p g6
##SBATCH -w n045,n046,n047,n048
#SBATCH --ntasks=6                # total number of nodes
#SBATCH --cpus-per-task=24         # total number of tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
## HPC ENVIRONMENT 
. /etc/profile.d/TMI.sh
##

mpirun ./pes > stdout.log
