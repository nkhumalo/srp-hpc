#!/bin/bash
#
echo "Begin"
#
mpirun -np 2 ../../../bin/LINUX64/nwchem 1h9t_prep_restrain.nw 2>&1 > 1h9t_prep_restrain.out
#
echo "Done Prepare"
#
qid=`qsub 1h9t_min_restrain.sh | sed 's/\.erato\.cs//' `
n=0
sleep 10
result=`qstat | grep $qid | grep -v C`
while [ ! -z "$result" ]
do
  n=$(($n+1))
  sleep 60
  result=`qstat | grep $qid | grep -v C`
  echo -n "$n..."
done
echo "$n"
cp 1h9t-restrain_md.qrs 1h9t-restrain_md.rst
#
echo "Done Minimize"
#
qid=`qsub 1h9t_eq_restrain.sh | sed 's/\.erato\.cs//' `
n=0
sleep 10
result=`qstat | grep $qid | grep -v C`
while [ ! -z "$result" ]
do
  n=$(($n+1))
  sleep 60
  result=`qstat | grep $qid | grep -v C`
  echo -n "$n..."
done
echo "$n"
#
echo "Done Equilibration"
#
mpirun -np 2 ../../../bin/LINUX64/nwchem 1h9t_eq_restrain_ana.nw 2>&1 > 1h9t_eq_restrain_ana.out
tar -zcf 1h9t_eq_restrain.tgz 1h9t_md.pdb 1h9t_eq_restrain.crd 1h9t-restrain_md.rst 1h9t-restrain.top
#
echo "Done Equilibration Analysis"
#
qid=`qsub 1h9t_md_restrain.sh | sed 's/\.erato\.cs//' `
n=0
sleep 10
result=`qstat | grep $qid | grep -v C`
while [ ! -z "$result" ]
do
  n=$(($n+1))
  sleep 120
  result=`qstat | grep $qid | grep -v C`
  echo -n "$n..."
done
echo "$n"
#
echo "Done Dynamics"
#
mpirun -np 2 ../../../bin/LINUX64/nwchem 1h9t_md_restrain_ana.nw 2>&1 > 1h9t_md_restrain_ana.out
tar -zcf 1h9t_md_restrain.tgz 1h9t_md.pdb 1h9t_md_restrain.crd
#
echo "Done Dynamics Analysis"
#
