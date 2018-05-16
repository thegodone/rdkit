#!/bin/bash
# pre-requiseted for boost & boost-pyrthon 1.65.1 build with c++11 option !!!
#brew install https://raw.githubusercontent.com/Homebrew/homebrew-core/7c9dc38d0749a863a696d87dda6d33e829799cb8/Formula/boost.rb
#brew install https://raw.githubusercontent.com/Homebrew/homebrew-core/7c9dc38d0749a863a696d87dda6d33e829799cb8/Formula/boost-python.rb
#ln -s /usr/local/Cellar/boost-python/1.65.1_1/lib/* /usr/local/Cellar/boost/1.65.1_1/lib/

PYROOT=/Library/Frameworks/Python.framework/Versions/2.7
SRCDIR=..
BLDDIR=build
#echo "Cleaning old directory"
rm -rf $BLDDIR
mkdir $BLDDIR
cd $BLDDIR
echo "cmake creating"
INSTALL=`pwd`/rdkit_build
BOOST_ROOT=/usr/local/Cellar/boost/1.65.1_1


CXXFLAGS="-std=c++11" \
cmake -DPYTHON_EXECUTABLE=$PYROOT/bin/python  \
          -DRDK_BUILD_AVALON_SUPPORT=ON \
          -DRDK_BUILD_INCHI_SUPPORT=ON \
          -DRDK_BUILD_CAIRO_SUPPORT=ON -DRDK_BUILD_THREADSAFE_SSS=ON -DRDK_TEST_MULTITHREADED=ON \
          -DPYTHON_LIBRARY=$PYROOT/lib/libpython2.7.dylib \
          -DPYTHON_INCLUDE_DIR=$PYROOT/include/python2.7 \
          -DPYTHON_NUMPY_INCLUDE_PATH=/usr/local/lib/python2.7/site-packages/numpy/core/include/ \
          -DCMAKE_BUILD_TYPE=Release \
          -DRDK_INSTALL_INTREE=OFF \
          -DCMAKE_INSTALL_PREFIX=$INSTALL \
          -DBOOST_ROOT=$BOOST_ROOT \
          -DBoost_NO_BOOST_CMAKE=FALSE \
          -DBoost_NO_SYSTEM_PATHS=on \
          -DBoost_USE_STATIC_LIBS=off \
          $SRCDIR && \
    make -j8 install && \
    RDBASE=`pwd`/$SRCDIR DYLD_FALLBACK_LIBRARY_PATH=`pwd`/rdkit_build/lib:$PYROOT/lib PYTHONPATH=`pwd`/rdkit_build/lib/python2.7/site-packages ctest -j8 \
    echo "To use this environment, do:" && \
    echo "export RDBASE="`pwd`/$SRCDIR && \
    echo "export DYLD_FALLBACK_LIBRARY_PATH="`pwd`/rdkit_build/lib:$PYROOT/lib && \
    echo "export PYTHONPATH=`pwd`/rdkit_build/lib/python2.7/site-packages"
# how to do it for matlab => you must do those two additional steps
# copy the rdkit_build/rdkit repertoire into /usr/local/lib/python2.7/site-packages/rdkit

#rm -rf usr/local/lib/python2.7/site-packages/rdkit/
#cp -rf  /Users/GVALMTGG/Github/rdkit/build/rdkit_build/lib/python2.7/site-packages/rdkit/ /usr/local/lib/python2.7/site-packages/rdkit

# symbolic link to be sure matlab can read the lib files from /usr/local/lib repertoire!

#rm -rf /usr/local/lib/libRDKit*
#ln -sf /Users/GVALMTGG/Github/rdkit/build/rdkit_build/lib/  /usr/local/lib

#export RDBASE=/Users/GVALMTGG/Github/rdkit/build/rdkit_build
#export PYROOT=/Library/Frameworks/Python.framework/Versions/2.7
#export DYLD_FALLBACK_LIBRARY_PATH=/Users/GVALMTGG/Github/rdkit/build/rdkit_build/lib:/Library/Frameworks/Python.framework/Versions/2.7/lib
#export PYTHONPATH=/Users/GVALMTGG/Github/rdkit/build/rdkit_build/lib/python2.7/site-packages
