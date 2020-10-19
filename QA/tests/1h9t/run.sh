#!/bin/bash
#
echo "Begin"
#
mpirun -np 2 ../../../bin/LINUX64/nwchem 1h9t_prep.nw 2>&1 > 1h9t_prep.out
#
echo "Done Prepare"
#
qid=`qsub 1h9t_min.sh | sed 's/\.erato\.cs//' `
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
cp 1h9t_md.qrs 1h9t_md.rst
#
echo "Done Minimize"
#
qid=`qsub 1h9t_eq.sh | sed 's/\.erato\.cs//' `
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
mpirun -np 2 ../../../bin/LINUX64/nwchem 1h9t_eq_ana.nw 2>&1 > 1h9t_eq_ana.out
tar -zcf 1h9t_eq.tgz 1h9t_md.pdb 1h9t_eq.crd
#
echo "Done Equilibration Analysis"
#
qid=`qsub 1h9t_md.sh | sed 's/\.erato\.cs//' `
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
mpirun -np 2 ../../../bin/LINUX64/nwchem 1h9t_md_ana.nw 2>&1 > 1h9t_md_ana.out
tar -zcf 1h9t_md.tgz 1h9t_md.pdb 1h9t_md.crd
#
echo "Done Dynamics Analysis"
#
