#!/bin/bash

#SBATCH --job-name=submit-seq-i.sh
#SBATCH -D .
#SBATCH --output=submit-seq-i.sh.o%j
#SBATCH --error=submit-seq-i.sh.e%j

USAGE="\n USAGE: ./submit-seq-i.sh prog\n
        prog        -> Program name\n"

if (test $# -lt 1 || test $# -gt 1)
then
        echo -e $USAGE
        exit 0
fi

make $1

export OMP_NUM_THREADS=1
export size=1073741824

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

export LD_PRELOAD=${EXTRAE_HOME}/lib/libomptrace.so
./$1 $size 
unset LD_PRELOAD

mpi2prv -f TRACE.mpits -o $1-1-${HOST}.prv -e $1 -paraver
rm -rf  TRACE.mpits set-0 >& /dev/null
