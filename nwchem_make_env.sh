# This file was auto-generated by /global/homes/l/ltang/nwchem-1/contrib/distro-tools/build_nwchem
export NWCHEM_TOP=/global/homes/l/ltang/nwchem-1
export NWCHEM_TARGET=LINUX64
export NWCHEM_MODULES="md"
export NWCHEM_MPIF_WRAP="/opt/cray/pe/craype/2.5.12/bin/ftn"
export NWCHEM_MPIC_WRAP="/opt/cray/pe/craype/2.5.12/bin/cc"
export NWCHEM_MPICXX_WRAP="/opt/cray/pe/craype/2.5.12/bin/CC"
export NWCHEM_LONG_PATHS=Y
export USE_NOFSCHECK=Y
export USE_MPI=y
export USE_MPIF=y
export USE_MPIF4=y
export MPI_INCLUDE=""
export MPI_LIB=""
export LIBMPI=""
export FC=ftn
export CC=cc
export CXX=CC
export ARMCI_NETWORK=MPI-TS
export MSG_COMMS=MPI
export SLURMOPT= 
export SLURM=
export BLASOPT="-lsci_intel_mpi -lsci_intel"
export BLAS_SIZE=4
export BLAS_LIB="-lsci_intel_mpi -lsci_intel `adios_config -l -f`"
export LAPACK_SIZE=4
export LAPACK_LIB="-lsci_intel_mpi -lsci_intel"
export SCALAPACK_SIZE=4
export SCALAPACK_LIB="-lsci_intel_mpi -lsci_intel"
export PYTHON_EXE=/usr/bin/python
export PYTHONVERSION=2.7
export USE_PYTHON64=yes
export PYTHONPATH=./:/global/homes/l/ltang/nwchem-1/contrib/python/
export PYTHONHOME=/usr
export PYTHONLIBTYPE=so
function renwc()
{
   make FC=$FC ; pushd $NWCHEM_TOP/src ; make FC=$FC link ; popd
}
export USE_64TO32=y
