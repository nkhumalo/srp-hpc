#!/bin/bash
#
# Optimize the parameters for the correlation functional
# using Monte Carlo.
#
#
declare -a splitter
#
export FCD00=0.198826
export FCD01=0.294438
export FCD11=0.247007
export FCD05=0.020225
export FCD51=-0.098842
export FCD55=-0.201236
#
export RCD00=0.884548
export RCD01=1.211432
export RCD11=0.255064
export RCD05=0.988291
export RCD51=0.535281
export RCD55=0.701049
#
export FCO00=-0.005884
export FCO01=0.063203
export FCO11=0.357171
export FCO05=0.019518
export FCO51=0.871725
export FCO55=1.051801
#
export RCO00=1.033050
export RCO01=0.850461
export RCO11=0.110139
export RCO05=0.985630
export RCO51=0.139729
export RCO55=0.918333
#
#
#
export DFCD00=0.1000
export DFCD01=0.1000
export DFCD11=0.1000
export DFCD05=0.1000
export DFCD51=0.1000
export DFCD55=0.1000
export DRCD00=0.100
export DRCD01=0.100
export DRCD11=0.100
export DRCD05=0.100
export DRCD51=0.100
export DRCD55=0.100
export DFCO00=0.1000
export DFCO01=0.1000
export DFCO11=0.1000
export DFCO05=0.1000
export DFCO51=0.1000
export DFCO55=0.1000
export DRCO00=0.100
export DRCO01=0.100
export DRCO11=0.100
export DRCO05=0.100
export DRCO51=0.100
export DRCO55=0.100
#
export factor_new=1.0
export factor_old=1.0
export shrink=0.5
count=0
misses=0
export REJECT="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
export ACCEPT="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
#
./nwchem-orbitals.sh "$FCD00 $FCD01 $FCD11 $FCD05 $FCD51 $FCD55 $RCD00 $RCD01 $RCD11 $RCD05 $RCD51 $RCD55  $FCO00 $FCO01 $FCO11 $FCO05 $FCO51 $FCO55 $RCO00 $RCO01 $RCO11 $RCO05 $RCO51 $RCO55"
parallel -a ./job_list.txt --colsep '\s+' -j ${PBS_NP} "./run_wf.sh {1} {2}"
./gather_results.sh
export AREA=`./compute_area.py`
#DEBUG
echo "HVD: AREA= $AREA" > junk.log
#DEBUG
echo "HVD: $count  $FCD00 $FCD01 $FCD11 $FCD05 $FCD51 $FCD55 $RCD00 $RCD01 $RCD11 $RCD05 $RCD51 $RCD55  $FCO00 $FCO01 $FCO11 $FCO05 $FCO51 $FCO55 $RCO00 $RCO01 $RCO11 $RCO05 $RCO51 $RCO55  $AREA"
echo "HVD: $count  $factor_old   $FCD00 $FCD01 $FCD11 $FCD05 $FCD51 $FCD55 $RCD00 $RCD01 $RCD11 $RCD05 $RCD51 $RCD55  $FCO00 $FCO01 $FCO11 $FCO05 $FCO51 $FCO55 $RCO00 $RCO01 $RCO11 $RCO05 $RCO51 $RCO55  $AREA" >  results_table.dat
echo "HVD: $count  $factor_old   $FCD00 $FCD01 $FCD11 $FCD05 $FCD51 $FCD55 $RCD00 $RCD01 $RCD11 $RCD05 $RCD51 $RCD55  $FCO00 $FCO01 $FCO11 $FCO05 $FCO51 $FCO55 $RCO00 $RCO01 $RCO11 $RCO05 $RCO51 $RCO55  $AREA" >  results_count.dat
#DEBUG
exit
#DEBUG
while true;
do
  export NEWSTEP=`./random.x $DFCD00 $DFCD01 $DFCD11 $DFCD05 $DFCD51 $DFCD55 $DRCD00 $DRCD01 $DRCD11 $DRCD05 $DRCD51 $DRCD55 $DFCO00 $DFCO01 $DFCO11 $DFCO05 $DFCO51 $DFCO55 $DRCO00 $DRCO01 $DRCO11 $DRCO05 $DRCO51 $DRCO55`
  export NEWVAR=`./add.x $FCD00 $FCD01 $FCD11 $FCD05 $FCD51 $FCD55 $RCD00 $RCD01 $RCD11 $RCD05 $RCD51 $RCD55  $FCO00 $FCO01 $FCO11 $FCO05 $FCO51 $FCO55 $RCO00 $RCO01 $RCO11 $RCO05 $RCO51 $RCO55 $NEWSTEP`
  ./nwchem-orbitals.sh "$NEWVAR"
  parallel -a ./job_list.txt --colsep '\s+' -j ${PBS_NP} "./run_wf.sh {1} {2}"
  ./gather_results.sh
  export AREA_T=`./compute_area.py`
  check=`echo "${AREA_T} < ${AREA}" | bc -l`
  count=`echo "$count + 1" | bc -l`
  echo "HVD p: $count $factor_old   $NEWVAR   $AREA_T" >> results_count.dat
  if [ $check -eq 1 ];
  then
    misses=0
    splitter=($NEWVAR)
    export AREA=$AREA_T
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
    export ACCEPT=`./max.x $ACCEPT $NEWSTEP`
    export NEWRNG=`./reshape.x $DFCD00 $DFCD01 $DFCD11 $DFCD05 $DFCD51 $DFCD55 $DRCD00 $DRCD01 $DRCD11 $DRCD05 $DRCD51 $DRCD55 $DFCO00 $DFCO01 $DFCO11 $DFCO05 $DFCO51 $DFCO55 $DRCO00 $DRCO01 $DRCO11 $DRCO05 $DRCO51 $DRCO55 $REJECT $ACCEPT`
    export REJECT="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
    export ACCEPT="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
    splitter=($NEWRNG)
    export DFCD00=${splitter[0]}
    export DFCD01=${splitter[1]}
    export DFCD11=${splitter[2]}
    export DFCD05=${splitter[3]}
    export DFCD51=${splitter[4]}
    export DFCD55=${splitter[5]}
    export DRCD00=${splitter[6]}
    export DRCD01=${splitter[7]}
    export DRCD11=${splitter[8]}
    export DRCD05=${splitter[9]}
    export DRCD51=${splitter[10]}
    export DRCD55=${splitter[11]}
    export DFCO00=${splitter[12]}
    export DFCO01=${splitter[13]}
    export DFCO11=${splitter[14]}
    export DFCO05=${splitter[15]}
    export DFCO51=${splitter[16]}
    export DFCO55=${splitter[17]}
    export DRCO00=${splitter[18]}
    export DRCO01=${splitter[19]}
    export DRCO11=${splitter[20]}
    export DRCO05=${splitter[21]}
    export DRCO51=${splitter[22]}
    export DRCO55=${splitter[23]}
    echo "HVD: $count  $NEWRNG"                      >> results_range.dat
    echo "HVD: $count  $NEWVAR  $AREA"
    echo "HVD: $count  $factor_old   $NEWVAR  $AREA" >> results_table.dat
  else
    export NEWVAR=`./subtract.x $FCD00 $FCD01 $FCD11 $FCD05 $FCD51 $FCD55 $RCD00 $RCD01 $RCD11 $RCD05 $RCD51 $RCD55  $FCO00 $FCO01 $FCO11 $FCO05 $FCO51 $FCO55 $RCO00 $RCO01 $RCO11 $RCO05 $RCO51 $RCO55 $NEWSTEP`
    ./nwchem-orbitals.sh "$NEWVAR"
    parallel -a ./job_list.txt --colsep '\s+' -j ${PBS_NP} "./run_wf.sh {1} {2}"
    ./gather_results.sh
    export AREA_T=`./compute_area.py`
    check=`echo "${AREA_T} < ${AREA}" | bc -l`
    echo "HVD m: $count   $factor_old   $NEWVAR   $AREA_T" >> results_count.dat
    if [ $check -eq 1 ];
    then
      misses=0
      splitter=($NEWVAR)
      export AREA=$AREA_T
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
      export ACCEPT=`./max.x $ACCEPT $NEWSTEP`
      export NEWRNG=`./reshape.x $DFCD00 $DFCD01 $DFCD11 $DFCD05 $DFCD51 $DFCD55 $DRCD00 $DRCD01 $DRCD11 $DRCD05 $DRCD51 $DRCD55 $DFCO00 $DFCO01 $DFCO11 $DFCO05 $DFCO51 $DFCO55 $DRCO00 $DRCO01 $DRCO11 $DRCO05 $DRCO51 $DRCO55 $REJECT $ACCEPT`
      export REJECT="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
      export ACCEPT="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
      splitter=($NEWRNG)
      export DFCD00=${splitter[0]}
      export DFCD01=${splitter[1]}
      export DFCD11=${splitter[2]}
      export DFCD05=${splitter[3]}
      export DFCD51=${splitter[4]}
      export DFCD55=${splitter[5]}
      export DRCD00=${splitter[6]}
      export DRCD01=${splitter[7]}
      export DRCD11=${splitter[8]}
      export DRCD05=${splitter[9]}
      export DRCD51=${splitter[10]}
      export DRCD55=${splitter[11]}
      export DFCO00=${splitter[12]}
      export DFCO01=${splitter[13]}
      export DFCO11=${splitter[14]}
      export DFCO05=${splitter[15]}
      export DFCO51=${splitter[16]}
      export DFCO55=${splitter[17]}
      export DRCO00=${splitter[18]}
      export DRCO01=${splitter[19]}
      export DRCO11=${splitter[20]}
      export DRCO05=${splitter[21]}
      export DRCO51=${splitter[22]}
      export DRCO55=${splitter[23]}
      echo "HVD: $count  $NEWRNG"                      >> results_range.dat
      echo "HVD: $count  $factor_old   $NEWVAR  $AREA"
      echo "HVD: $count  $factor_old   $NEWVAR  $AREA" >> results_table.dat
    else
      misses=`echo "$misses + 1" | bc -l`
      export REJECT=`./max.x $REJECT $NEWSTEP`
    fi
  fi
  check_misses=`echo "$misses >= 10" | bc -l`
  if [ $check_misses -eq 1 ];
  then
    #
    # Many consecutive misses so shrink the Monte Carlo box
    #
    misses=0
    export NEWRNG=`./scale.x $DFCD00 $DFCD01 $DFCD11 $DFCD05 $DFCD51 $DFCD55 $DRCD00 $DRCD01 $DRCD11 $DRCD05 $DRCD51 $DRCD55  $DFCO00 $DFCO01 $DFCO11 $DFCO05 $DFCO51 $DFCO55 $DRCO00 $DRCO01 $DRCO11 $DRCO05 $DRCO51 $DRCO55 $shrink`
    splitter=($NEWRNG)
    export DFCD00=${splitter[0]}
    export DFCD01=${splitter[1]}
    export DFCD11=${splitter[2]}
    export DFCD05=${splitter[3]}
    export DFCD51=${splitter[4]}
    export DFCD55=${splitter[5]}
    export DRCD00=${splitter[6]}
    export DRCD01=${splitter[7]}
    export DRCD11=${splitter[8]}
    export DRCD05=${splitter[9]}
    export DRCD51=${splitter[10]}
    export DRCD55=${splitter[11]}
    export DFCO00=${splitter[12]}
    export DFCO01=${splitter[13]}
    export DFCO11=${splitter[14]}
    export DFCO05=${splitter[15]}
    export DFCO51=${splitter[16]}
    export DFCO55=${splitter[17]}
    export DRCO00=${splitter[18]}
    export DRCO01=${splitter[19]}
    export DRCO11=${splitter[20]}
    export DRCO05=${splitter[21]}
    export DRCO51=${splitter[22]}
    export DRCO55=${splitter[23]}
    export factor_old=$factor_new
    export factor_new=`echo "$factor_new * $shrink" | bc -l`
  fi
  check_factor=`echo "$factor_old < 0.00001" | bc -l`
  if [ $check_factor -eq 1 ];
  then
    exit
  fi
done
