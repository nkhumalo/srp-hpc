#!/bin/bash
export NWCHEM_EXE=/home/hvandam/nwchem-1-wfn1/bin/LINUX64/nwchem
#if [ -f results_wf.dat ]
#then
#  rm results_wf.dat
#fi
if [ -f job_list.txt ]
then
  rm job_list.txt
fi
if [ -f job_list_6_31g.txt ]
then
  rm job_list_6_31g.txt
fi
if [ -f job_list_sto_3g.txt ]
then
  rm job_list_sto_3g.txt
fi
if [ -f results_hf_6_31g.dat ]
then
  rm results_hf_6_31g.dat
fi
if [ -f results_hf_sto_3g.dat ]
then
  rm results_hf_sto_3g.dat
fi
if [ -f results_mc_6_31g.dat ]
then
  rm results_mc_6_31g.dat
fi
if [ -f results_mc_sto_3g.dat ]
then
  rm results_mc_sto_3g.dat
fi

R=1.138692
#R=10.988692
check=`echo "$R < 111.3" | bc -l`
R12=`echo "$R / 2.000000" | bc -l`
while [ $check -eq 1 ]
do
  #rm mcscf_dat.*
  cat << EOF > input-mc-6-31g.nw
echo
start mcscf_dat
geometry units au noautoz
  H 0 0 -$R12
  H 0 0  $R12
  symmetry cs
end
basis
  H library 6-31g 
end
mcscf
  active       4
  actelec      2
  multiplicity 1
end
task mcscf
EOF
  #$NWCHEM_EXE input-mc-6-31g.nw 2>&1 > nwchem-mc-6-31g-${R}.out
  energy_mc_6_31g=`grep "Total MCSCF energy" nwchem-mc-6-31g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm mcscf_dat.*
  cat << EOF > input-mc-sto-3g.nw
echo
start mcscf_dat
geometry units au noautoz
  H 0 0 -$R12
  H 0 0  $R12
  symmetry cs
end
basis
  H library sto-3g
end
mcscf
  active       2
  actelec      2
  multiplicity 1
end
task mcscf
EOF
  #$NWCHEM_EXE input-mc-sto-3g.nw 2>&1 > nwchem-mc-sto-3g-${R}.out
  energy_mc_sto_3g=`grep "Total MCSCF energy" nwchem-mc-sto-3g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm hf_dat.*
  cat << EOF > input-hf-6-31g.nw
echo
start hf_dat
geometry units au noautoz
  H 0 0 -$R12
  H 0 0  $R12
  symmetry cs
end
basis
  H library 6-31g 
end
task scf energy
EOF
  #$NWCHEM_EXE input-hf-6-31g.nw 2>&1 > nwchem-hf-6-31g-${R}.out
  energy_hf_6_31g=`grep "Total SCF energy" nwchem-hf-6-31g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm hf_dat.*
  cat << EOF > input-hf-sto-3g.nw
echo
start hf_dat
geometry units au noautoz
  H 0 0 -$R12
  H 0 0  $R12
  symmetry cs
end
basis
  H library sto-3g
end
task scf energy
EOF
  #$NWCHEM_EXE input-hf-sto-3g.nw 2>&1 > nwchem-hf-sto-3g-${R}.out
  energy_hf_sto_3g=`grep "Total SCF energy" nwchem-hf-sto-3g-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm wfn1_dat.*
  echo "input-wf-6-31g-${R}.nw nwchem-wf-6-31g-${R}.out ${R}" >> job_list.txt
  echo "input-wf-6-31g-${R}.nw nwchem-wf-6-31g-${R}.out ${R}" >> job_list_6_31g.txt
  cat << EOF > input-wf-6-31g-${R}.nw
echo
start wfn1_dat
geometry units au noautoz
  H 0 0 -$R12
  H 0 0  $R12
  symmetry cs
end
basis
  H library 6-31g 
end
set wfn1:input_vectors "atomic"
set wfn1:solver "monte-carlo"
set wfn1:print_error F
set wfn1:t_initial 0.1
#set wfn1:corr_expr "wfn1_mx"
set wfn1:corr_expr "wfn1c6"
set wfn1:fac_ab2a  $1
set wfn1:fac_ab2b  $2
set wfn1:fac_ab2c  $3
set wfn1:fac_ab4a  $4
set wfn1:fac_ab4b  $5
set wfn1:fac_ab4c  $6
set wfn1:t_bath 0.0
set wfn1:maxit  20000
task wfn1 energy
EOF
  echo "input-wf-sto-3g-${R}.nw nwchem-wf-sto-3g-${R}.out ${R}" >> job_list.txt
  echo "input-wf-sto-3g-${R}.nw nwchem-wf-sto-3g-${R}.out ${R}" >> job_list_sto_3g.txt
  cat << EOF > input-wf-sto-3g-${R}.nw
echo
start wfn1_dat
geometry units au noautoz
  H 0 0 -$R12
  H 0 0  $R12
  symmetry cs
end
basis
  H library sto-3g
end
set wfn1:input_vectors "atomic"
set wfn1:solver "monte-carlo"
set wfn1:print_error F
set wfn1:t_initial 0.1
#set wfn1:corr_expr "wfn1_mx"
set wfn1:corr_expr "wfn1c6"
set wfn1:fac_ab2a  $1
set wfn1:fac_ab2b  $2
set wfn1:fac_ab2c  $3
set wfn1:fac_ab4a  $4
set wfn1:fac_ab4b  $5
set wfn1:fac_ab4c  $6
set wfn1:t_bath 0.0
set wfn1:maxit  20000
task wfn1 energy
EOF
  echo $R $energy_hf_6_31g  >> results_hf_6_31g.dat
  echo $R $energy_hf_sto_3g >> results_hf_sto_3g.dat
  echo $R $energy_mc_6_31g  >> results_mc_6_31g.dat
  echo $R $energy_mc_sto_3g >> results_mc_sto_3g.dat
  R=`echo "$R + 0.25" | bc -l`
  check=`echo "$R < 9.0" | bc -l`
  R12=`echo "$R / 2.000000" | bc -l`
done
cat job_list_6_31g.txt job_list_sto_3g.txt > job_list.txt
