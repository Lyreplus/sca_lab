#!/bin/bash

JOB_SCRIPT="submit-mpi.sh"
SIZE=12000
ARRAY_A=( 1 2 4 8 16 32 64 )
ARRAY_B=( 1 2 3 )

TEMP_JOB_FILE="temp_job.sh"
cp $JOB_SCRIPT $TEMP_JOB_FILE

for A in "${ARRAY_A[@]}"; do
        sed -i "s/^#SBATCH --output=.*/#SBATCH --output=output_${A}/" $TEMP_JOB_FILE
        # sed -i "s/^export size=.*/export size=${A}/" $TEMP_JOB_FILE
        sed -i "s/^#SBATCH --ntasks=.*/#SBATCH --ntasks=1/" $TEMP_JOB_FILE
        sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=${A}/" $TEMP_JOB_FILE
        # sed -i "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=${NTASKS}/" $TEMP_JOB_FILE
        sbatch $TEMP_JOB_FILE $A
        sleep 15
done

rm $TEMP_JOB_FILE

FILE="out.txt"

echo "All jobs have been submitted."
for A in "${ARRAY_A[@]}"; do
        echo "${A} = $(cat "output_${A}")" >> $FILE
done

cat $FILE