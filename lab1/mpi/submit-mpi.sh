#!/bin/bash

#SBATCH --job-name=submit-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-mpi.sh.o%j
#SBATCH --error=submit-mpi.sh.e%j
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1

PROGRAM=pi_mpi
size=1073741824
#size=$1

srun --mpi=pmi2 --cpu-bind=socket ./$PROGRAM $size

