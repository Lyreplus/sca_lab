#!/bin/bash
export SOFT_ROOT=/Soft/PAR
export KMP_AFFINITY=scatter

# FFTW3 library
export FFTW3_HOME=$SOFT_ROOT/fftw3
export LD_LIBRARY_PATH=$FFTW3_HOME/lib

# Extrae
export EXTRAE_HOME=$SOFT_ROOT/extrae/current
#export EXTRAE_HOME=$SOFT_ROOT/extrae/previous
export PATH=$EXTRAE_HOME/bin/:$PATH
export LD_LIBRARY_PATH=/Soft/libunwind/current/lib:$EXTRAE_HOME/lib:$LD_LIBRARY_PATH
export EXTRAE_ON=TRUE
export EXTRAE_CONFIG_FILE=$SOFT_ROOT/extrae.xml
#export EXTRAE_CONFIG_FILE=/scratch/nas/1/par0/djg/modelfactors/lab1/extrae.xml

# Paraver
export PARAVER_HOME=$SOFT_ROOT/paraver/current
export PATH=$PARAVER_HOME/bin:$PATH

# Dimemas for Tareador simulations
export TAREADOR_HOME=$SOFT_ROOT/TAREADOR
export PATH=$TAREADOR_HOME/install/trf2trf/bin:$PATH
export PATH=$TAREADOR_HOME/src/scripts/dimemas_simulations:$PATH
#export PATH=$TAREADOR_HOME/install/dimemas-5.2.3-linux-x86_64/bin:$PATH
export DIMEMAS_HOME=$SOFT_ROOT/dimemas/current
export PATH=$DIMEMAS_HOME/bin:$PATH

# BSC's basic analysis tools
export BASICANALYSIS_HOME=$SOFT_ROOT/basicanalysis/current
export PATH=$BASICANALYSIS_HOME:$PATH
#export CLUSTERING_HOME=$HOME/clusteringsuite
#export PATH=$CLUSTERING_HOME/bin:$PATH

# Tareador GuiED
#export TAREADOR_GUI_ROOT=/scratch/nas/1/pap0/TAREADOR/gui
export TAREADOR_GUI_ROOT=$TAREADOR_HOME/src/TAREADOR_GUI_LLVM
export TAREADOR_ICONS=$TAREADOR_GUI_ROOT/xdot/icons
export PATH=$TAREADOR_GUI_ROOT/backend:$TAREADOR_GUI_ROOT/xdot:$PATH
export PYTHONPATH=$TAREADOR_GUI_ROOT:$PYTHONPATH

# Tareador LLVM
export TAREADOR_LLVM_HOME=$TAREADOR_HOME/install/llvm-clang
export TAREADOR_INSTALL_DIRECTORY=$TAREADOR_LLVM_HOME
export PATH=$TAREADOR_LLVM_HOME/bin:$PATH

# Tareador compilation and execution environment
alias tareador='schroot -p -c tareador'

# To get rid of problems ...
export LC_ALL="en_US.utf-8"
alias evince="LD_LIBRARY_PATH= /usr/bin/evince"

# Intel compiler and tools
#source /Soft/intel_parallel_studio_xe/parallel_studio_xe_2020/psxevars.sh
#source /Soft/PAR/psxevars.sh
export PATH=/Soft/gcc/10.2.0/bin:$PATH
export FI_PROVIDER=sockets
