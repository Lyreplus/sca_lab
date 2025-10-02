#!/bin/bash

JOB_SCRIPT="submit-omp.sh"
ARRAY_A=( 625 1250 2500 5000 10000 20000 40000 )
TEMP_JOB_FILE="temp_job.sh"

cp $JOB_SCRIPT $TEMP_JOB_FILE

# Array per salvare gli ID dei job
JOB_IDS=()

for A in "${ARRAY_A[@]}"; do
        sed -i "s/^#SBATCH --output=.*/#SBATCH --output=output_${A}/" $TEMP_JOB_FILE
        sed -i "s/^export size=.*/export size=${A}/" $TEMP_JOB_FILE

        # Sottometti il job e salva l'ID
        JOB_ID=$(sbatch $TEMP_JOB_FILE BackSubs_omp 32 | awk '{print $4}')
        JOB_IDS+=($JOB_ID)
        echo "Submitted job $JOB_ID for size ${A}"

        sleep 10
done

rm $TEMP_JOB_FILE

echo "All jobs have been submitted."
echo "Waiting for all jobs to complete..."

# Aspetta che tutti i job siano completati
for JOB_ID in "${JOB_IDS[@]}"; do
    while squeue -j $JOB_ID 2>/dev/null | grep -q $JOB_ID; do
            sleep 30
    done
    echo "Job $JOB_ID completed"
done

echo "All jobs completed!"
echo ""
echo "Results:"

for A in "${ARRAY_A[@]}"; do
    if [ -f "output_${A}" ]; then
            echo "${A} = $(cat "output_${A}" | grep "elapsed")"
    else
            echo "${A} = Output file not found"
    fi
done

