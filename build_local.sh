#! /bin/bash
TARGET=$1
export INSTALL_PREFIX=$(echo $CONDA_PREFIX)
rm -rf build
mkdir build
cd build
CC=/usr/bin/gcc CXX=/usr/bin/g++ cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX ..
make $TARGET
cd ..
