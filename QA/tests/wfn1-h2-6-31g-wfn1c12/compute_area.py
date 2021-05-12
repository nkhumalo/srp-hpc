#!/bin/env python3
#
# Compute the area of differences between the MCSCF and WFN1
# results. 
#
import sys
import math
#
f_mcscf=open("results_mc_6_31g.dat")
l_mcscf=f_mcscf.readlines()
f_mcscf.close()
f_mcscf=open("results_mc_sto_3g.dat")
l_mcscf.extend(f_mcscf.readlines())
f_mcscf.close()
#
f_wfn1=open("results_wf_6_31g.dat")
l_wfn1=f_wfn1.readlines()
f_wfn1.close()
f_wfn1=open("results_wf_sto_3g.dat")
l_wfn1.extend(f_wfn1.readlines())
f_wfn1.close()
#
area=0.0
len_mcscf=len(l_mcscf)
len_wfn1=len(l_wfn1)
if len_mcscf != len_wfn1:
  print("Mismatching file lengths:",len_mcscf,len_wfn1)
  sys.exit()
for ii in range(0,len_wfn1):
  line_mcscf = l_mcscf[ii]
  line_wfn1  = l_wfn1[ii]
  p_mcscf    = line_mcscf.split()
  p_wfn1     = line_wfn1.split()
  x_mcscf    = float(p_mcscf[0])
  x_wfn1     = float(p_wfn1[0])
  e_mcscf    = float(p_mcscf[1])
  e_wfn1     = float(p_wfn1[1])
  if x_mcscf != x_wfn1:
    print("x-coordinates out of synch",ii,x_mcscf,x_wfn1)
    sys.exit()
  #area += abs(e_mcscf-e_wfn1)*abs(e_mcscf-e_wfn1)
  area = max(area,abs(e_mcscf-e_wfn1))
#area = math.sqrt(area)
print(area)
