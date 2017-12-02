#- env: == default ==
os=`uname`
export NWCHEM_TOP=$TRAVIS_BUILD_DIR
export USE_MPI=y
export USE_64TO32=y
 if [[ "$os" == "Darwin" ]]; then 
export BLASOPT="-L/usr/local/opt/openblas/lib -lopenblas"
export SCALAPACK="-L/usr/local/lib  -lscalapack -lopenblas"
export NWCHEM_TARGET=MACX64 
export PATH=/usr/local/opt/open-mpi/bin/:$PATH 
export DYLD_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$DYLD_LIBRARY_PATH
fi
if [[ "$os" == "Linux" ]]; then 
   export BLASOPT="-L$TRAVIS_BUILD_DIR/lib -lopenblas"
   export SCALAPACK="-L$TRAVIS_BUILD_DIR/lib  -lscalapack -lopenblas"
   export NWCHEM_TARGET=LINUX64 
   export LD_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$LD_LIBRARY_PATH
fi
export BLAS_SIZE=4
export SCALAPACK_SIZE=4
export USE_PYTHONCONFIG=y
export PYTHONVERSION=2.7
export PYTHONHOME=/usr
