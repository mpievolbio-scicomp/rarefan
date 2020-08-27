#! /bin/sh
# source me!
if [ ! -d $CONDA_PREFIX ]
then
    echo "Activate the conda environment first."
fi

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
