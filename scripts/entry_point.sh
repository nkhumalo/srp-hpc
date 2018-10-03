#!/bin/bash

# TODO: Improve entry point to work with various cases and folder mount options
if [ ! -f ~/.nwchemrc ]; then 
  cd QA; ./domknwchemrc; cd /nwchem
fi
cd QA/tests/ethanol; 
mpirun -np 4 ../../../bin/LINUX64/nwchem ethanol_md.nw
