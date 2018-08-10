#!/bin/env python
import os
import glob
import subprocess
#
# This script runs the dynamics and analysis "workflow". Both the dynamics and
# the analysis are executed concurrently. To get this to work we first launch 
# NWChem to run the MD simulation in the nt-6-6-5-md directory. Then we loop
# until the trajectory file appears. At that point we can start the analysis run.
#
os.chdir("./nt-6-6-5-md")
if os.path.isfile("nt-ice_md.trj"):
    os.remove("nt-ice_md.trj")
#
# Run the "prepare" step and cleanup trace and profile files afterwards
# (the performance of this step is not important)
#
subprocess.call("../../../../bin/LINUX64/nwchem nt-6-6-5-a-ice-1x1x1-opt-oh-h-prep.nw | tee nt-6-6-5-a-ice-1x1x1-opt-oh-h-prep.out",shell=True)
files = glob.glob("tautrace.*")+glob.glob("profile.*")+glob.glob("events.*")
for file in files:
    os.remove(file)
#
# Start the "dynamics" step in the background and wait for the trajectory file to 
# appear
#
file = open("nt-6-6-5-a-ice-1x1x1-opt-oh-h-md.out","w")
subprocess.call("mpirun -np 4 ../../../../bin/LINUX64/nwchem nt-6-6-5-a-ice-1x1x1-opt-oh-h-md.nw 2>&1 > nt-6-6-5-a-ice-1x1x1-opt-oh-h-md.out &",shell=True)
while (not os.path.isfile("nt-ice_md.trj")):
    pass
#
os.chdir("../nt-6-6-5-ana")
#
# Start the "analysis"
#
subprocess.call("../../../../bin/LINUX64/nwchem nt-6-6-5-a-ice-1x1x1-opt-oh-h-ana.nw | tee nt-6-6-5-a-ice-1x1x1-opt-oh-h-ana.out",shell=True)
