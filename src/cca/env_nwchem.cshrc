# NWChem settings
setenv NWCHEM_TOP /home/windus/testnwchem/nwchem
setenv NWCHEM_TARGET LINUX
setenv NWCHEM_MODULES "nwdft gradients"
#setenv USE_SUBGROUPS 2
setenv USE_SHARED y

# GA settings
setenv USE_MPI y
setenv TARGET LINUX

# MPI settings
setenv MPI_PATH /msrc/apps/mpich-1.2.6/gcc/ch_shmem
setenv MPI_LIB $MPI_PATH/lib
setenv MPI_INCLUDE $MPI_PATH/include
setenv LIBMPI -lmpich
setenv PATH $MPI_PATH/bin:$PATH

#MCMD Settings 
# NOTE: Makesure all the env settings from INSTALL file in the following
#      directory (/home/vidhya/CCA/mcmd-manoj/MCMD-NWCHEM), is specified here
#
#setenv GA_PATH /home/windus/CCA/mcmd-paper/nwchem-shared-sg/src/tools
setenv GA_PATH /home/windus/testnwchem/nwchem/src/tools
#setenv NWCHEM_CCA_ROOT /home/windus/CCA/mcmd-paper/sacomp/nwchem

setenv PATH /home/windus/CCA/cca-tools-0.6.0_rc1/local/bin:$PATH
