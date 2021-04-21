#!/bin/bash
#
# Execute NWChem on a single input file
#
# The command has two arguments
# - $1 is the input file
# - $2 is the output file
#
# The scripts proceeds as follows:
# - create the scratch directory
# - change into the scratch directory
# - run NWChem writing the output to the output file
# - cat the occupation numbers to the end of the output file
# - clean the scratch directory up
#
export PATH=`readlink -f ~/bin`:$PATH
export MK_SCRATCH=`readlink -f ~/bin/test_make_scratch_serial.sh`
export MK_CLEAN=`readlink -f ~/bin/test_make_clean_serial.sh`
export NWCHEM_EXE=/home/hvandam/nwchem-1-wfn1/bin/LINUX64/nwchem

. ${MK_SCRATCH}

${NWCHEM_EXE} ${CMS_WORKDIR}/$1 2>&1 > ${CMS_WORKDIR}/$2
cat fort.80 >> ${CMS_WORKDIR}/$2

. ${MK_CLEAN}
