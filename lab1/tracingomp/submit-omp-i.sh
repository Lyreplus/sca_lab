#!/bin/bash

#SBATCH --job-name=submit-omp-i.sh
#SBATCH -D .
#SBATCH --output=submit-omp-i.sh.o%j
#SBATCH --error=submit-omp-i.sh.e%j

USAGE="\n USAGE: ./submit-omp-i.sh prog numthreads \n
        prog        -> Program name\n
        numthreads  -> Number of threads in parallel execution\n"

if (test $# -lt 2 || test $# -gt 2)
then
        echo -e $USAGE
        exit 0
fi

make $1

export OMP_NUM_THREADS=$2
export size=1073741824

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

export LD_PRELOAD=${EXTRAE_HOME}/lib/libomptrace.so
./$1 $size
unset LD_PRELOAD

mpi2prv -f TRACE.mpits -o $1-$2-${HOST}.prv -e $1 -paraver
rm -rf  TRACE.mpits set-0 >& /dev/null
