#!/bin/bash

. /etc/profile


cd /Codar/nwchem-1/QA/tests/ethanol

cp /Codar/nwchem-1/contrib/codar_integration/md_compress_adaptive_adios_dataspaces.py .
cp /Codar/nwchem-1/contrib/codar_integration/dataspaces.conf .

sed -i 's/coord 0/coord 1/' ethanol_md.nw
sed -i 's/scoor 0/scoor 1/' ethanol_md.nw


mpirun -n 1 dataspaces_server -s 1 -c 2 &

sleep 10

mpirun -n 1 ../../../bin/LINUX64/nwchem ethanol_md.nw &

sleep 5

mpirun -n 1 python3 ./md_compress_adaptive_adios_dataspaces.py &

wait

