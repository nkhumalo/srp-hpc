#!/bin/bash
export NWCHEM_EXE=../../../bin/LINUX64/nwchem
#if [ -f results_wf.dat ]
#then
#  rm results_wf.dat
#fi
if [ -f job_list.txt ]
then
  rm job_list.txt
fi
if [ -f results_hf.dat ]
then
  rm results_hf.dat
fi
if [ -f results_mc.dat ]
then
  rm results_mc.dat
fi

R=1.138692
#R=10.988692
check=`echo "$R < 111.3" | bc -l`
R12=`echo "$R / 2.000000" | bc -l`
while [ $check -eq 1 ]
do
  #rm mcscf_dat.*
  cat << EOF > input-mc.nw
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
  #$NWCHEM_EXE input-mc.nw 2>&1 > nwchem-mc-${R}.out
  energy_mc=`grep "Total MCSCF energy" nwchem-mc-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm hf_dat.*
  cat << EOF > input-hf.nw
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
  #$NWCHEM_EXE input-hf.nw 2>&1 > nwchem-hf-${R}.out
  energy_hf=`grep "Total SCF energy" nwchem-hf-${R}.out | tail -n 1 | awk '{ print $5 }'`
  #rm wfn1_dat.*
  echo "input-wf-${R}.nw nwchem-wf-${R}.out ${R}" >> job_list.txt
  cat << EOF > input-wf-${R}.nw
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
set wfn1:corr_expr "wfn1c2"
set wfn1:fac_ab2   $1
set wfn1:fac_ab4a  $2
set wfn1:fac_ab4b  $3
set wfn1:fac_ab4c  $4
set wfn1:t_bath 0.0
set wfn1:maxit  20000
task wfn1 energy
EOF
  #energy_wf=`grep "Total WFN1 energy" nwchem.out | tail -n 1 | awk '{ print $5 }'`
  #occs=`cat fort.80`
  #echo $R $energy_wf $occs >> results_wf.dat
  echo $R $energy_hf       >> results_hf.dat
  echo $R $energy_mc       >> results_mc.dat
  R=`echo "$R + 0.25" | bc -l`
  check=`echo "$R < 9.0" | bc -l`
  R12=`echo "$R / 2.000000" | bc -l`
done
