#!/bin/bash
#set -x
########################################################################
#   Batch Job Parameters
########################################################################
#
#### Grid Engine Job Wrapper ####
#
#$ -N 1h9t_md_restrain
#$ -pe mpi2 1
#$ -q wp08
#$ -j y
#$ -m e
#$ -v LD_LIBRARY_PATH,OMP_NUM_THREADS
#
#### Torque/PBS Job Parameters ####
#
#PBS -N 1h9t_md_restrain
#PBS -l nodes=1:ppn=16
#
########################################################################
#   Batch Environment Initialization Commands
########################################################################
#
#### Grid Engine Initialization Commands ####
#
if [ -f $TMPDIR/sge_init.sh ]; then
    source $TMPDIR/sge_init.sh
fi
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
export OMP_NUM_THREADS=1
export GMON_OUT_PREFIX="gmon_md"
#echo "LD_LIBRARY_PATH="$LD_LIBRARY_PATH
echo "PWD"
pwd
echo "LS"
ls -l /home/hvandam/nwchem-1/QA/tests/1h9t/1h9t_md_restrain.nw
echo "NWCHEM"
mpiexec -n 16 -x LD_LIBRARY_PATH /home/hvandam/nwchem-1/bin/LINUX64/nwchem /home/hvandam/nwchem-1/QA/tests/1h9t/1h9t_md_restrain.nw 2>&1 > 1h9t_md_restrain.out
