#!/bin/bash
declare -a files
#
if [ -f results_h2_wf_6_31g.dat ];
then
  rm results_h2_wf_6_31g.dat
fi
while IFS= read -r line; do
  files=($line)
  output=${files[1]}
  xcoord=${files[2]}
  energy=`grep "Total WFN1 energy" ${output} | tail -n 1 | awk '{ print $5 }'`
  occs=`tail -n 1 ${output}`
  echo $xcoord $energy $occs >> results_h2_wf_6_31g.dat
done < "./job_list_h2_6_31g.txt"
#
if [ -f results_h2_wf_sto_3g.dat ];
then
  rm results_h2_wf_sto_3g.dat
fi
while IFS= read -r line; do
  files=($line)
  output=${files[1]}
  xcoord=${files[2]}
  energy=`grep "Total WFN1 energy" ${output} | tail -n 1 | awk '{ print $5 }'`
  occs=`tail -n 1 ${output}`
  echo $xcoord $energy $occs >> results_h2_wf_sto_3g.dat
done < "./job_list_h2_sto_3g.txt"
#
if [ -f results_lih_wf_sto_3g.dat ];
then
  rm results_lih_wf_sto_3g.dat
fi
while IFS= read -r line; do
  files=($line)
  output=${files[1]}
  xcoord=${files[2]}
  energy=`grep "Total WFN1 energy" ${output} | tail -n 1 | awk '{ print $5 }'`
  occs=`tail -n 1 ${output}`
  echo $xcoord $energy $occs >> results_lih_wf_sto_3g.dat
done < "./job_list_lih_sto_3g.txt"
