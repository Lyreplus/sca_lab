#!/bin/bash

#SBATCH --job-name=submit-seq.sh
#SBATCH -D .
#SBATCH --output=submit-seq.sh.o%j
#SBATCH --error=submit-seq.sh.e%j


USAGE="\n USAGE: submit-seq.sh PROG size\n
        PROG   -> program name\n
        size   -> size of the problem\n"

if (test $# -lt 2 || test $# -gt 2)
then
	echo -e $USAGE
    exit 0
fi

make $1

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

if (test "${HOST}" = "boada-11")
then
    echo "Use sbatch to execute this script"
    exit 0
elif (test "${HOST}" = "boada-12")
then
    echo "Use sbatch to execute this script"
    exit 0
elif (test "${HOST}" = "boada-13")
then
    echo "Use sbatch to execute this script"
    exit 0
fi

/usr/bin/time -o time-$1-${HOST} ./$1 $2
