#!/bin/bash
export NWCHEM_EXE=../../../bin/LINUX64/nwchem
if [ -f results_wf.dat ]
then
  rm results_wf.dat
fi
if [ -f results_hf.dat ]
then
  rm results_hf.dat
fi
if [ -f results_mc.dat ]
then
  rm results_mc.dat
fi

R=2.18778447
#R=10.988692
check=`echo "$R < 111.3" | bc -l`
R12=`echo "$R / 1.000000" | bc -l`
while [ $check -eq 1 ]
do
  rm mcscf_dat.*
  cat << EOF > input-mc.nw
echo
start mcscf_dat
geometry units au noautoz
  H  0 0 -$R12
  Be 0 0  0
  H  0 0  $R12
  symmetry cs
end
basis
  H  library 6-31g
  Be library 6-31g
end
mcscf
  active       13
  actelec      6
  multiplicity 1
end
task mcscf
EOF
  $NWCHEM_EXE input-mc.nw 2>&1 > nwchem-mc-${R}.out
  energy_mc=`grep "Total MCSCF energy" nwchem-mc-${R}.out | tail -n 1 | awk '{ print $5 }'`
  rm hf_dat.*
  cat << EOF > input-hf.nw
echo
start hf_dat
geometry units au noautoz
  H  0 0 -$R12
  Be 0 0  0
  H  0 0  $R12
  symmetry cs
end
basis
  H  library 6-31g
  Be library 6-31g
end
task scf energy
EOF
  $NWCHEM_EXE input-hf.nw 2>&1 > nwchem-hf-${R}.out
  energy_hf=`grep "Total SCF energy" nwchem-hf-${R}.out | tail -n 1 | awk '{ print $5 }'`
  rm wfn1_dat.*
  cat << EOF > input-wf.nw
echo
start wfn1_dat
geometry units au noautoz
  H  0 0 -$R12
  Be 0 0  0
  H  0 0  $R12
  symmetry cs
end
basis
  H  library 6-31g
  Be library 6-31g
end
set wfn1:input_vectors "atomic"
set wfn1:solver "monte-carlo"
set wfn1:print_error F
set wfn1:t_initial 0.1
#set wfn1:corr_expr "wfn1_mx"
set wfn1:corr_expr "wfn1c"
set wfn1:t_bath 0.0
set wfn1:maxit  20000
task wfn1 energy
EOF
  $NWCHEM_EXE input-wf.nw 2>&1 > nwchem.out
  cp nwchem.out nwchem-wf-${R}.out
  energy_wf=`grep "Total WFN1 energy" nwchem.out | tail -n 1 | awk '{ print $5 }'`
  occs=`cat fort.80`
  orbe_hf=`sort -u fort.70`
  orbe_mc=`cat     fort.71`
  #orbe_na=`cat     fort.75`
  #orbe_ca=`cat     fort.76`
  #orbe_hfna=`cat     fort.77`
  #orbe_hfca=`cat     fort.78`
  #orbe_hfpa=`cat     fort.79`
  #orbe_junk=`cat     fort.85`
  #orbe_crna=`cat     fort.89`
  #orbe_crca=`cat     fort.90`
  #orbe_dea=`cat     fort.110`
  #orbe_deb=`cat     fort.111`
  echo $R $orbe_hf      >> orb_hf_energies.dat
  echo $R $orbe_mc      >> orb_mc_energies.dat
  #echo $R $orbe_na      >> orb_na_energies.dat
  #echo $R $orbe_ca      >> orb_ca_energies.dat
  #echo $R $orbe_hfna    >> orb_hfna_energies.dat
  #echo $R $orbe_hfca    >> orb_hfca_energies.dat
  #echo $R $orbe_hfpa    >> orb_hfpa_energies.dat
  #echo $R $orbe_junk    >> orb_junk_energies.dat
  #echo $R $orbe_crna    >> orb_crna_energies.dat
  #echo $R $orbe_crca    >> orb_crca_energies.dat
  #echo $R $orbe_dea     >> orb_dea_energies.dat
  #echo $R $orbe_deb     >> orb_deb_energies.dat
  echo $R $energy_wf $occs >> results_wf.dat
  echo $R $energy_hf       >> results_hf.dat
  echo $R $energy_mc       >> results_mc.dat
  echo $R $energy_wf $occs
  check=`echo "$R > 11.2" | bc -l`
  if [ $check -eq 1 ]; then
    R=`echo "$R + 1000000.0" | bc -l`
  else
    R=`echo "$R + 0.1" | bc -l`
  fi
  check=`echo "$R < 1000011.3" | bc -l`
  check=`echo "$R < 11.3" | bc -l`
  #DEBUG
  #check=`echo "$R < 1.6" | bc -l`
  #DEBUG
  R12=`echo "$R / 1.000000" | bc -l`
done
