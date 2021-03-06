#!/bin/bash
export NWCHEM_EXE=../../../bin/LINUX64/nwchem
$NWCHEM_EXE input-hf.nw 2>&1 > result_hf.out
$NWCHEM_EXE input-mc.nw 2>&1 > result_mc.out
$NWCHEM_EXE input-wf.nw 2>&1 > result_wf.out
