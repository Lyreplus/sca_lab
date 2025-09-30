#!/bin/bash

JOB_SCRIPT="submit-mpi.sh"
NTASKS_ARRAY=(  1 2 3 4 5 6 7 8 9 10 \
                11 12 13 14 15 16 17 18 19 20 \
                21 22 23 24 25 26 27 28 29 30 
                31 32 33 34 35 36 37 38 39 40 \
                41 42 43 44 45 46 47 48 49 50 \
                51 52 53 54 55 56 57 58 59 60 \
                61 62 63 64 )
NNODES_ARRAY=( 1 2 3 )

for NNODES in "${NNODES_ARRAY[@]}"; do
        for NTASKS in "${NTASKS_ARRAY[@]}"; do
                TEMP_JOB_FILE="temp_job.sh"
                cp $JOB_SCRIPT $TEMP_JOB_FILE
                sed -i "s/^#SBATCH --output=.*/#SBATCH --output=output_mpi_${NTASKS}_${NNODES}/" $TEMP_JOB_FILE
                sed -i "s/^#SBATCH --ntasks=.*/#SBATCH --ntasks=${NNODES}/" $TEMP_JOB_FILE
                sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=${NNODES}/" $TEMP_JOB_FILE
                sed -i "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=${NTASKS}/" $TEMP_JOB_FILE
                sbatch $TEMP_JOB_FILE
                sleep 10
                rm $TEMP_JOB_FILE
        done
done

echo "All jobs have been submitted."
for NNODES in "${NNODES_ARRAY[@]}"; do
        for NTASKS in "${NTASKS_ARRAY[@]}"; do
                echo "${NNODES} ${NTASKS} = $(cat "output_mpi_${NTASKS}_${NNODES}")"
        done
done