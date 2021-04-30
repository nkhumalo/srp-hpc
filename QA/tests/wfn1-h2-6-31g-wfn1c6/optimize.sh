#!/bin/bash
#
# Optimize the parameters for the correlation functional
# using Monte Carlo.
#
declare -a splitter
export FAB2A=0.3062
export FAB2B=0.3352
export FAB2C=0.0133
export FAB4A=0.7109
export FAB4B=0.3865
export FAB4C=0.5318
export RAB2A=0.01
export RAB2B=0.01
export RAB2C=0.01
export RAB4A=0.01
export RAB4B=0.01
export RAB4C=0.01
count=0
./nwchem-orbitals.sh $FAB2A $FAB2B $FAB2C $FAB4A $FAB4B $FAB4C
parallel -a ./job_list.txt --colsep '\s+' -j ${PBS_NP} "./run_wf.sh {1} {2}"
./gather_results.sh
export AREA=`./compute_area.py`
echo "HVD: $FAB2A $FAB2B $FAB2C $FAB4A $FAB4B $FAB4C $AREA"
echo "HVD: $count $FAB2A $FAB2B $FAB2C $FAB4A $FAB4B $FAB4C $AREA" >> results_table.dat
echo "HVD: $count $FAB2A $FAB2B $FAB2C $FAB4A $FAB4B $FAB4C $AREA" >  results_count.dat
while true;
do
  export NEWVAR=`./random.x $FAB2A $FAB2B $FAB2C $FAB4A $FAB4B $FAB4C $RAB2A $RAB2B $RAB2C $RAB4A $RAB4B $RAB4C`
  ./nwchem-orbitals.sh $NEWVAR
  parallel -a ./job_list.txt --colsep '\s+' -j ${PBS_NP} "./run_wf.sh {1} {2}"
  ./gather_results.sh
  export AREA_T=`./compute_area.py`
  check=`echo "${AREA_T} < ${AREA}" | bc -l`
  count=`echo "$count + 1" | bc -l`
  echo "HVD: $count $NEWVAR $AREA_T" >> results_count.dat
  if [ $check -eq 1 ];
  then
    splitter=($NEWVAR)
    export AREA=$AREA_T
    export FAB2A=${splitter[0]}
    export FAB2B=${splitter[1]}
    export FAB2C=${splitter[2]}
    export FAB4A=${splitter[3]}
    export FAB4B=${splitter[4]}
    export FAB4C=${splitter[5]}
    echo "HVD: $FAB2A $FAB2B $FAB2C $FAB4A $FAB4B $FAB4C $AREA"
    echo "HVD: $count $FAB2A $FAB2B $FAB2C $FAB4A $FAB4B $FAB4C $AREA" >> results_table.dat
  fi
done
