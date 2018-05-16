#!/bin/bash
PYROOT=/Users/GVALMTGG/miniconda2
SRCDIR=../
BLDDIR=build_condapython27_release
#echo "Cleaning old directory"
cd $BLDDIR
echo "cmake creating"
INSTALL=`pwd`/rdkit_build
BOOST_ROOT=$PYROOT

#RDBASE=`pwd`/$SRCDIR DYLD_FALLBACK_LIBRARY_PATH=`pwd`/rdkit_build/lib:$PYROOT/lib PYTHONPATH=`pwd`/rdkit_build/lib/python2.7/site-packages ctest -j4 \
echo "To use this environment, do:" && \
echo "export RDBASE="`pwd`/$SRCDIR && \
echo "export DYLD_FALLBACK_LIBRARY_PATH="`pwd`/rdkit_build/lib:$PYROOT/lib && \
echo "export PYTHONPAH=`pwd`/rdkit_build/lib/python2.7/site-packages"