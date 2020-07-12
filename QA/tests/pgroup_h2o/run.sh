#!/bin/bash
mpirun -np 12 ../../../bin/LINUX64/nwchem pgroup_h2o.nw | tee pgroup_h2o.out
