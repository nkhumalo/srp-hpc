#!/bin/bash
#PBS -N nt-6-6-5-md
#PBS -l nodes=1:ppn=16
#
mpirun -np 16 /home/hvandam/nwchem-1/bin/LINUX64/nwchem /home/hvandam/nwchem-1/QA/tests/nt-6-6-5-10h2o/nt-6-6-5-md.nw 2>&1 > /home/hvandam/nwchem-1/QA/tests/nt-6-6-5-10h2o/nt-6-6-5-md.out
