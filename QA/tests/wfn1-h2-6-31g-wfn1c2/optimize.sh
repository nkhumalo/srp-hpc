#!/bin/bash
#
# Optimize the parameters for the correlation functional
# using Monte Carlo.
#
declare -a splitter
export FAB2=0.17419223781349386
export FAB4A=1.1148553279846742
export FAB4B=1.2408066078452848
export FAB4C=0.60550732663181817
export RAB2=0.025
export RAB4A=0.025
export RAB4B=0.05
export RAB4C=0.05
./nwchem-orbitals.sh $FAB2 $FAB4A $FAB4B $FAB4C
parallel -a ./job_list.txt --colsep '\s+' -j ${PBS_NP} "./run_wf.sh {1} {2}"
./gather_results.sh
export AREA=`./compute_area.py`
echo "HVD: $FAB2 $FAB4A $FAB4B $FAB4C $AREA"
#DEBUG
exit
#DEBUG
while true;
do
  export NEWVAR=`./random.x $FAB2 $FAB4A $FAB4B $FAB4C $RAB2 $RAB4A $RAB4B $RAB4C`
  ./nwchem-orbitals.sh $NEWVAR
  parallel -a ./job_list.txt --colsep '\s+' -j ${PBS_NP} "./run_wf.sh {1} {2}"
  ./gather_results.sh
  export AREA_T=`./compute_area.py`
  check=`echo "${AREA_T} < ${AREA}" | bc -l`
  if [ $check -eq 1 ];
  then
    splitter=($NEWVAR)
    export AREA=$AREA_T
    export FAB2=${splitter[0]}
    export FAB4A=${splitter[1]}
    export FAB4B=${splitter[2]}
    export FAB4C=${splitter[3]}
    echo "HVD: $FAB2 $FAB4A $FAB4B $FAB4C $AREA"
  fi
done
