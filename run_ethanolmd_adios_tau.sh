#!/bin/bash

export TAU_PROFILE=1
export TAU_TRACE=1

. /etc/profile

cd QA/tests/ethanol

sed -i 's/coord 0/coord 1/' ethanol_md.nw
sed -i 's/scoor 0/scoor 1/' ethanol_md.nw

export PYTHONPATH=${PYTHONPATH}:/opt/adios2/lib/python3.5/site-packages
export PYTHONPATH=${PYTHONPATH}:/MDTrAnal/lib
mpirun --allow-run-as-root -n 10 python3 /MDTrAnal/src/md_compress_adaptive_adios2.py &
sleep 10

mpirun --allow-run-as-root -n 2 ../../../bin/LINUX64/nwchem ethanol_md.nw

#bpls -l nwchem_xyz.bp

