#!/bin/bash

JOB_SCRIPT="submit-seq.sh"
ARRAY_A=( 5000 10000 20000 40000 80000 160000 )

TEMP_JOB_FILE="temp_job.sh"
cp $JOB_SCRIPT $TEMP_JOB_FILE

for A in "${ARRAY_A[@]}"; do
        sed -i "s/^#SBATCH --output=.*/#SBATCH --output=output_${A}/" $TEMP_JOB_FILE
        sed -i "s/^export size=.*/export size=${A}/" $TEMP_JOB_FILE
        # sed -i "s/^#SBATCH --ntasks=.*/#SBATCH --ntasks=${NNODES}/" $TEMP_JOB_FILE
        # sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=${NNODES}/" $TEMP_JOB_FILE
        # sed -i "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=${NTASKS}/" $TEMP_JOB_FILE
        sbatch $TEMP_JOB_FILE
        sleep 10
done

rm $TEMP_JOB_FILE

echo "All jobs have been submitted."
for A in "${ARRAY_A[@]}"; do
        echo "${A} = $(cat "output_${A}")"
done