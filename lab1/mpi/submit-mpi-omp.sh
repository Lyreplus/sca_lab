#!/bin/bash

#SBATCH --job-name=submit-mpi-omp.sh
#SBATCH -D .
#SBATCH --output=submit-mpi-omp.sh.o%j
#SBATCH --error=submit-mpi-omp.sh.e%j
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=12

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

PROGRAM=pi_mpi_omp
size=1073741824
#size=$1

srun --mpi=pmi2 --cpu-bind=socket ./$PROGRAM $size
