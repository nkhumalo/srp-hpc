#!/bin/bash
#
# TAU_CXX_ABSOLUTE: The absolute path name of the C++ compiler. This is required
# to avoid an infinite recursion arising from calling this script.
#
# TAU_MAKEFILE: The absolute path of a modified Tau makefile. This makefile needs
# to have been modified to use the TAU_*_ABSOLUTE compilers instead of looking for
# the compiler in the PATH.
#
if [ ${#USE_TAU} -eq 0 ]; then
  ${TAU_CXX_ABSOLUTE} -g $@
else
  # TAU_MAKEFILE has been hacked to ensure absolute pathnames for the compilers
  #export TAU_MAKEFILE=/global/homes/v/vandam/bin/Makefile.tau-intel-mpi-pdt
  if [ ${USE_TAU} == "source" ]; then
    # export F90=/opt/cray/pe/craype/2.5.12/bin/ftn
    tau_f90.sh -g $@
  elif [ ${USE_TAU} == "compiler" ]; then
    # export F90=/opt/cray/pe/craype/2.5.12/bin/ftn
    tauf90 -g $@
  else
    echo "Unknown USE_TAU: " $USE_TAU
    exit 10
  fi
fi
