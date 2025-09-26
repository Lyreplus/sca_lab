#!/bin/bash

#SBATCH --job-name=submit-arch.sh
#SBATCH -D .
#SBATCH --output=submit-arch.sh.o%j
#SBATCH --error=submit-arch.sh.e%j

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

if (test "${HOST}" = "boada-6")
then
    echo "Use sbatch to execute this script"
    exit 0
elif (test "${HOST}" = "boada-7")
then
    echo "Use sbatch to execute this script"
    exit 0
elif (test "${HOST}" = "boada-8")
then
    echo "Use sbatch to execute this script"
    exit 0
fi

PROG=lscpu
$PROG > ${PROG}-${HOST}

PROG='lstopo'
$PROG > ${PROG}-${HOST}
PROGFIG='lstopo --of fig map.fig'
$PROGFIG
mv map.fig map-${HOST}.fig
