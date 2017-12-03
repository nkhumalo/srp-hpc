#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 if [[ "$NWCHEM_MODULES" == "tce" ]]; then
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 tce_n2 tce_ccsd_t_h2o tce_h2o_eomcc
   cd $NWCHEM_TOP/QA/testoutputs
   diff tce_h2o_eomcc.o*nw*
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 tce_ipccsd_f2
   diff tce_ipccsd_f2.o*nw*

 else
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 dft_siosi3
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs  2 h2o_opt dft_he2+ cosmo_h2o_dft tddft_h2o h2o2-response
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs  2 pspw 
 fi
