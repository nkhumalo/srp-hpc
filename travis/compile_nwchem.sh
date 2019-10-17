#!/bin/bash
echo "start compile"
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
ls -lrt $TRAVIS_BUILD_DIR|tail -3
os=`uname`
arch=`uname -m`
if [[ "$NWCHEM_MODULES" == "tce" ]]; then 
    export EACCSD=1
    export IPCCSD=1
fi
cd $TRAVIS_BUILD_DIR/src
if [[ "$arch" == "aarch64" ]]; then 
    if [[ "$NWCHEM_MODULES" == "tce" ]]; then 
	FOPT2="-O0 -fno-aggressive-loop-optimizations"
    else
	FOPT2="-O1 -fno-aggressive-loop-optimizations"
    fi
else    
    FOPT2="-O2 -fno-aggressive-loop-optimizations"
fi    
 if [[ "$os" == "Darwin" ]]; then 
   if [[ "$NWCHEM_MODULES" == "tce" ]]; then
     FOPT2="-O1 -fno-aggressive-loop-optimizations"
   fi
   if [[ ! -z "$USE_SIMINT" ]] ; then 
       FOPT2="-O0 -fno-aggressive-loop-optimizations"
       SIMINT_BUILD_TYPE=Debug
   fi
     ../travis/sleep_loop.sh make  FDEBUG="-O0 -g" FOPTIMIZE="$FOPT2" -j3
     cd $TRAVIS_BUILD_DIR/src/64to32blas 
     make
     cd $TRAVIS_BUILD_DIR/src
     ../contrib/getmem.nwchem 1000
     otool -L ../bin/MACX64/nwchem
#     printenv DYLD_LIBRARY_PATH
#     ls -lrt $DYLD_LIBRARY_PATH
#      tail -120 make.log
 elif [[ "$os" == "Linux" ]]; then
     if [[ "$arch" == "aarch64" ]]; then 
	 export MAKEFLAGS=-j8
     else    
	 export MAKEFLAGS=-j3
     fi
     ../travis/sleep_loop.sh make  FDEBUG="-O0 -g" FOPTIMIZE="$FOPT2" 
     cd $TRAVIS_BUILD_DIR/src/64to32blas 
     make
     cd $TRAVIS_BUILD_DIR/src
     $TRAVIS_BUILD_DIR/contrib/getmem.nwchem 1000
 fi
 #caching
 mkdir -p $TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET $TRAVIS_BUILD_DIR/.cachedir/files
 cp $TRAVIS_BUILD_DIR/bin/$NWCHEM_TARGET/nwchem  $NWCHEM_EXECUTABLE
 echo === ls binaries cache ===
 ls -lrt $TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/ || true
 echo =========================
 rsync -av $TRAVIS_BUILD_DIR/src/basis/libraries  $TRAVIS_BUILD_DIR/.cachedir/files/.
