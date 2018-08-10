#!/bin/bash
#set -x
#
#### Torque/PBS Job Parameters ####
#
#PBS -N wfn1s
#PBS -l nodes=16
#
########################################################################
#   Batch Environment Initialization Commands
########################################################################
#
#### Torque/PBS Initialization Commands ####
#
if [ ! -z $PBS_JOBID ]; then
    cd $PBS_O_WORKDIR
fi
#
########################################################################
#   Job Execution Commands
########################################################################
#
export OMP_NUM_THREADS=16

make --jobs=16

