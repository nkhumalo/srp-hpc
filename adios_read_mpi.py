#!/usr/bin/env python
"""
Example:

$ python ./adios_read_mpi.py [METHOD] [params]
"""

import adios_mpi as ad
import getopt, sys
import os
import numpy as np

method = "BP"
init = "verbose=3;check_read_status=0"

if len(sys.argv) > 1:
    method = sys.argv[1]

if len(sys.argv) > 2:
    init = sys.argv[2]

print(">>> Method:", method, init)
ad.read_init(method, parameters=init)

f = ad.file("/Codar/nwchem-1/QA/tests/ethanol/tau-metrics.bp", method, is_stream=True, timeout_sec = -1.0)
print(f)

i = 0
while True:
    print(">>> step:", i)
    v = f.var['counter_values'][...]
    print('counter_values=', np.sum(v))
    if ('program_name 0' in f.attr):
        a = f.attr['program_name 0'][...]
        print('program_name 0=', a)

    if (f.advance(timeout_sec=-1.0) < 0):
        break
    i += 1

f.close()

ad.read_finalize(method)

print(">>> Done.")
