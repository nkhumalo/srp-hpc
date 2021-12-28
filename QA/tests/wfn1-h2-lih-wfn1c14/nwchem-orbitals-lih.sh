#!/bin/bash
#
# We have many command line arguments (24 or more). This does not work with
# bash. So all command line arguments are packed into 1. We first need to
# split that using an array, before we can do something sensible.
#
declare -a splitter
export NWCHEM_EXE=/home/hvandam/nwchem-1-wfn1/bin/LINUX64/nwchem-new
#if [ -f results_wf.dat ]
#then
#  rm results_wf.dat
#fi
#if [ -f job_list.txt ]
#then
#  rm job_list.txt
#fi
#if [ -f job_list_lih_6_31g.txt ]
#then
#  rm job_list_lih_6_31g.txt
#fi
if [ -f job_list_lih_sto_3g.txt ]
then
  rm job_list_lih_sto_3g.txt
fi
#if [ -f results_lih_hf_6_31g.dat ]
#then
#  rm results_lih_hf_6_31g.dat
#fi
if [ -f results_lih_hf_sto_3g.dat ]
then
  rm results_lih_hf_sto_3g.dat
fi
#if [ -f results_lih_mc_6_31g.dat ]
#then
#  rm results_lih_mc_6_31g.dat
#fi
if [ -f results_lih_mc_sto_3g.dat ]
then
  rm results_lih_mc_sto_3g.dat
fi
#echo "$@" > called_nwchem_orbitals.txt
splitter=($1)
export FCD00=${splitter[0]}
export FCD01=${splitter[1]}
export FCD11=${splitter[2]}
export FCD05=${splitter[3]}
export FCD51=${splitter[4]}
export FCD55=${splitter[5]}
export RCD00=${splitter[6]}
export RCD01=${splitter[7]}
export RCD11=${splitter[8]}
export RCD05=${splitter[9]}
export RCD51=${splitter[10]}
export RCD55=${splitter[11]}
export FCO00=${splitter[12]}
export FCO01=${splitter[13]}
export FCO11=${splitter[14]}
export FCO05=${splitter[15]}
export FCO51=${splitter[16]}
export FCO55=${splitter[17]}
export RCO00=${splitter[18]}
export RCO01=${splitter[19]}
export RCO11=${splitter[20]}
export RCO05=${splitter[21]}
export RCO51=${splitter[22]}
export RCO55=${splitter[23]}
#
R=2.65904199
#R=10.988692
check=`echo "$R < 111.3" | bc -l`
R12=`echo "$R / 2.000000" | bc -l`
while [ $check -eq 1 ]
do
  #rm mcscf_dat.*
  cat << EOF > input-lih-mc-6-31g-${R}.nw
echo
start mcscf_dat
geometry units au noautoz
  Li 0 0 -$R12
  H  0 0  $R12
  symmetry c1
end
basis
  Li library 6-31g 
  H  library 6-31g 
end
scf
  vectors input atomic output hf
end
task scf
mcscf
  active       11
  actelec      4
  multiplicity 1
  vectors input hf output mcscf
end
task mcscf
EOF
  #mpirun -np 16 $NWCHEM_EXE input-lih-mc-6-31g-${R}.nw 2>&1 > nwchem-lih-mc-6-31g-${R}.out
  #energy_lih_mc_6_31g=`grep "Total MCSCF energy" nwchem-lih-mc-6-31g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm mcscf_dat.*
  cat << EOF > input-lih-mc-sto-3g-${R}.nw
echo
start mcscf_dat
geometry units au noautoz
  Li 0 0 -$R12
  H  0 0  $R12
  symmetry c1
end
basis
  Li library sto-3g
  H  library sto-3g
end
scf
  vectors input atomic output hf
end
task scf
mcscf
  active       6
  actelec      4
  multiplicity 1
  vectors input hf output mcscf
end
task mcscf
EOF
  #mpirun -np 1 $NWCHEM_EXE input-lih-mc-sto-3g-${R}.nw 2>&1 > nwchem-lih-mc-sto-3g-${R}.out
  energy_lih_mc_sto_3g=`grep "Total MCSCF energy" nwchem-lih-mc-sto-3g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm hf_dat.*
  cat << EOF > input-lih-hf-6-31g-${R}.nw
echo
start hf_dat
geometry units au noautoz
  Li 0 0 -$R12
  H  0 0  $R12
  symmetry c1
end
basis
  Li library 6-31g 
  H  library 6-31g 
end
task scf energy
EOF
  #$NWCHEM_EXE input-lih-hf-6-31g-${R}.nw 2>&1 > nwchem-lih-hf-6-31g-${R}.out
  #energy_lih_hf_6_31g=`grep "Total SCF energy" nwchem-lih-hf-6-31g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm hf_dat.*
  cat << EOF > input-lih-hf-sto-3g-${R}.nw
echo
start hf_dat
geometry units au noautoz
  Li 0 0 -$R12
  H  0 0  $R12
  symmetry c1
end
basis
  Li library sto-3g
  H  library sto-3g
end
task scf energy
EOF
  #$NWCHEM_EXE input-lih-hf-sto-3g-${R}.nw 2>&1 > nwchem-lih-hf-sto-3g-${R}.out
  energy_lih_hf_sto_3g=`grep "Total SCF energy" nwchem-lih-hf-sto-3g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm wfn1_dat.*
  #echo "input-lih-wf-6-31g-${R}.nw nwchem-lih-wf-6-31g-${R}.out ${R}" >> job_list_lih.txt
  #echo "input-lih-wf-6-31g-${R}.nw nwchem-lih-wf-6-31g-${R}.out ${R}" >> job_list_lih_6_31g.txt
  cat << EOF > input-lih-wf-6-31g-${R}.nw
echo
start wfn1_dat
geometry units au noautoz
  Li 0 0 -$R12
  H  0 0  $R12
  symmetry c1
end
basis
  Li library 6-31g 
  H  library 6-31g 
end
set wfn1:input_vectors "atomic"
set wfn1:solver "monte-carlo"
set wfn1:print_error F
set wfn1:t_initial 0.1
#set wfn1:corr_expr "wfn1_mx"
set wfn1:corr_expr "wfn1c14"
#
set wfn1:fcd00  $FCD00
set wfn1:fcd01  $FCD01
set wfn1:fcd11  $FCD11
set wfn1:fcd05  $FCD05
set wfn1:fcd51  $FCD51
set wfn1:fcd55  $FCD55
#
set wfn1:rcd00  $RCD00
set wfn1:rcd01  $RCD01
set wfn1:rcd11  $RCD11
set wfn1:rcd05  $RCD05
set wfn1:rcd51  $RCD51
set wfn1:rcd55  $RCD55
#
set wfn1:fco00  $FCO00
set wfn1:fco01  $FCO01
set wfn1:fco11  $FCO11
set wfn1:fco05  $FCO05
set wfn1:fco51  $FCO51
set wfn1:fco55  $FCO55
#
set wfn1:rco00  $RCO00
set wfn1:rco01  $RCO01
set wfn1:rco11  $RCO11
set wfn1:rco05  $RCO05
set wfn1:rco51  $RCO51
set wfn1:rco55  $RCO55
#
set wfn1:t_bath 0.0
set wfn1:maxit  20000
task wfn1 energy
EOF
  #echo "input-lih-wf-sto-3g-${R}.nw nwchem-lih-wf-sto-3g-${R}.out ${R}" >> job_list.txt
  echo "input-lih-wf-sto-3g-${R}.nw nwchem-lih-wf-sto-3g-${R}.out ${R}" >> job_list_lih_sto_3g.txt
  cat << EOF > input-lih-wf-sto-3g-${R}.nw
echo
start wfn1_dat
geometry units au noautoz
  Li 0 0 -$R12
  H  0 0  $R12
  symmetry c1
end
basis
  Li library sto-3g
  H  library sto-3g
end
set wfn1:input_vectors "atomic"
set wfn1:solver "monte-carlo"
set wfn1:print_error F
set wfn1:t_initial 0.1
#set wfn1:corr_expr "wfn1_mx"
set wfn1:corr_expr "wfn1c14"
#
set wfn1:fcd00  $FCD00
set wfn1:fcd01  $FCD01
set wfn1:fcd11  $FCD11
set wfn1:fcd05  $FCD05
set wfn1:fcd51  $FCD51
set wfn1:fcd55  $FCD55
#
set wfn1:rcd00  $RCD00
set wfn1:rcd01  $RCD01
set wfn1:rcd11  $RCD11
set wfn1:rcd05  $RCD05
set wfn1:rcd51  $RCD51
set wfn1:rcd55  $RCD55
#
set wfn1:fco00  $FCO00
set wfn1:fco01  $FCO01
set wfn1:fco11  $FCO11
set wfn1:fco05  $FCO05
set wfn1:fco51  $FCO51
set wfn1:fco55  $FCO55
#
set wfn1:rco00  $RCO00
set wfn1:rco01  $RCO01
set wfn1:rco11  $RCO11
set wfn1:rco05  $RCO05
set wfn1:rco51  $RCO51
set wfn1:rco55  $RCO55
#
set wfn1:t_bath 0.0
set wfn1:maxit  20000
task wfn1 energy
EOF
  #echo $R $energy_lih_hf_6_31g  >> results_lih_hf_6_31g.dat
  echo $R $energy_lih_hf_sto_3g >> results_lih_hf_sto_3g.dat
  #echo $R $energy_lih_mc_6_31g  >> results_lih_mc_6_31g.dat
  echo $R $energy_lih_mc_sto_3g >> results_lih_mc_sto_3g.dat
  R=`echo "$R + 0.25" | bc -l`
  check=`echo "$R < 12.0" | bc -l`
  R12=`echo "$R / 2.000000" | bc -l`
done
#cat job_list_6_31g.txt job_list_sto_3g.txt > job_list.txt
cat job_list_lih_sto_3g.txt job_list_h2_6_31g.txt job_list_h2_sto_3g.txt > job_list.txt
