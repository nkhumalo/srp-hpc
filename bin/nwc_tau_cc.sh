#!/bin/bash
#
# TAU_CC_ABSOLUTE: The absolute path name of the C compiler. This is required
# to avoid an infinite recursion arising from calling this script.
#
# TAU_MAKEFILE: The absolute path of a modified Tau makefile. This makefile needs
# to have been modified to use the TAU_*_ABSOLUTE compilers instead of looking for
# the compiler in the PATH.
#
# However, hacking TAU_MAKEFILE does not work as someone must have decided to be
# "clever" and take the basename of the absolute name provided for the compiler.
# Hence to terminate the recursion the USE_BRAIN variable is introduced. If this 
# variable is set we call the compiler by it's absolute pathname, no matter what.
#
if [ ${#USE_BRAIN} -eq 0 ]; then
  export USE_BRAIN="yes"
  if [ ${#USE_TAU} -eq 0 ]; then
    ${TAU_CC_ABSOLUTE} -g $@
  else
    # TAU_MAKEFILE has been hacked to ensure absolute pathnames for the compilers
    #export TAU_MAKEFILE=/global/homes/v/vandam/bin/Makefile.tau-intel-mpi-pdt
    if [ ${USE_TAU} == "source" ]; then
      # export F90=/opt/cray/pe/craype/2.5.12/bin/ftn
      tau_cc.sh -g $@
    elif [ ${USE_TAU} == "compiler" ]; then
      # export F90=/opt/cray/pe/craype/2.5.12/bin/ftn
      taucc -g $@
    else
      echo "Unknown USE_TAU: " $USE_TAU
      exit 10
    fi
  fi
else
  ${TAU_CC_ABSOLUTE} -g $@
fi
