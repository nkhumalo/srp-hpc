#- env: == default ==
os=`uname`
arch=`uname -m`
#export NWCHEM_TOP=$TRAVIS_BUILD_DIR
export USE_MPI=y
if [[ "$os" == "Darwin" ]]; then 
   if [[ "$NWCHEM_MODULES" = "tce" ]]; then
     export USE_INTERNALBLAS=y
   else
     export USE_64TO32="y"
     export BLASOPT="-L/usr/local/opt/openblas/lib -lopenblas"
     export LAPACK_LIB="-L/usr/local/opt/openblas/lib -lopenblas"
     if [[ "$MPI_IMPL" == "openmpi" ]]; then
       export SCALAPACK="-L/usr/local/lib -lscalapack -lopenblas"
     fi
   fi
  export NWCHEM_TARGET=MACX64 
  export DYLD_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$DYLD_LIBRARY_PATH
  if [[ "$MPI_IMPL" == "openmpi" ]]; then
    export PATH=/usr/local/opt/open-mpi/bin/:$PATH 
  fi
  if [[ "$MPI_IMPL" == "mpich" ]]; then 
    export PATH=/usr/local/opt/mpich/bin/:$PATH 
  fi
fi
if [[ "$os" == "Linux" ]]; then 
   export NWCHEM_TARGET=LINUX64 
   if [[ -z "$USE_SIMINT" ]] && [[ "$arch" != "aarch64" ]]; then 
     export BUILD_OPENBLAS="y"
     export BUILD_SCALAPACK="y"
   else
     export BLASOPT="-lopenblas"
     export LAPACK_LIB="-lopenblas"
     if [[ "$MPI_IMPL" == "mpich" ]]; then 
       export SCALAPACK="-lscalapack-mpich -lopenblas"
     elif [[ "$MPI_IMPL" == "openmpi" ]]; then
       export SCALAPACK="-lscalapack-openmpi -lopenblas"
     fi
   fi
     export USE_64TO32="y"
#   fi
fi
export OMP_NUM_THREADS=1
export USE_NOIO=1
if [[ "$USE_64TO32" == "y" ]]; then
  export BLAS_SIZE=4
  export SCALAPACK_SIZE=4
fi
#export USE_PYTHONCONFIG=y
#export PYTHONVERSION=2.7
#export PYTHONHOME=/usr
