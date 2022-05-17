#! /bin/bash
cat /proc/meminfo
conda init bash
source ~/.bashrc
conda activate repinpop
export INSTALL_PREFIX=$(echo $CONDA_PREFIX)
mkdir build
cd build
CC=/usr/bin/gcc CXX=/usr/bin/g++ cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX ..
cd ..
cd REPIN_ecology/REPIN_ecology 
gradle build
cd ../..
