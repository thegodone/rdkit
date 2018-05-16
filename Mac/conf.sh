#!/bin/bash
# pre-requiseted for boost & boost-pyrthon 1.65.1 build with c++11 option !!!
#brew install https://raw.githubusercontent.com/Homebrew/homebrew-core/7c9dc38d0749a863a696d87dda6d33e829799cb8/Formula/boost.rb
#brew install https://raw.githubusercontent.com/Homebrew/homebrew-core/7c9dc38d0749a863a696d87dda6d33e829799cb8/Formula/boost-python.rb
#ln -s /usr/local/Cellar/boost-python/1.65.1_1/lib/* /usr/local/Cellar/boost/1.65.1_1/lib/

PYROOT=/Library/Frameworks/Python.framework/Versions/2.7
SRCDIR=../
BLDDIR=build_python27_release
#echo "Cleaning old directory"

echo "cmake creating"
INSTALL=`pwd`/rdkit_build
BOOST_ROOT=/usr/local/Cellar/boost/1.65.1_1



A=RDBASE=`pwd`/$SRCDIR DYLD_FALLBACK_LIBRARY_PATH=`pwd`/rdkit_build/lib:$PYROOT/lib PYTHONPATH=`pwd`/rdkit_build/lib/python2.7/site-packages
echo A
echo "To use this environment, do:" && \
echo "export RDBASE="`pwd`/$SRCDIR && \
echo "export DYLD_FALLBACK_LIBRARY_PATH="`pwd`/rdkit_build/lib:$PYROOT/lib && \
echo "export PYTHONPATH=`pwd`/rdkit_build/lib/python2.7/site-packages"