# RepinPop

## Compilers and build system
Make sure, a C++ compiler and the `autoconf` utility are installed on your system.

## Create the conda environment

```
$> conda env create -n repinpop --file=environment.yml
```

This step will create a new conda environment named `repinpop` and install a number of dependencies from the conda repository. For
details, inspect the file `environment.yml` which lists the dependencies and
version requirements. 

Activate the new environment:

```
$> conda activate repinpop
```
Record the value of the `$CONDA_PREFIX` environment variable, e.g.

```
echo $CONDA_PREFIX > conda_prefix.txt
```


## Install 3rd party libraries through cmake.
Not all dependencies are available on the conda archives. `andi` [Efficient
Estimation of Evolutionary Distances](https://github.com/EvolBioInf/andi.git),
it's dependency `divsufsort` [A lightweight suffix-sorting library](https://github.com/y-256/libdivsufsort.git), and `clustDist` [Cluster Distances into Phylogenies](https://github.com/EvolBioInf/clustDist.git) are handled by a the script `CMakeLists.txt` to be consumed by the `cmake` utility. But first, we have to deactivate the conda environment. Make sure you have GSL
installed in your system or add the GSL libraries and header file locations to your environment.
```
conda deactivate
```
```
$> mkdir build
$> cd build
$> cmake -DCMAKE_INSTALL_PREFIX=`cat conda_prefix.txt` ..
```
The last line instructs cmake to setup the `$CONDA_PREFIX` as the installation prefix for the third party libraries to be installed.
```
 $> make
```

This will download the required source codes for all three dependencies, build,
and install the executables into the `conda` environment created in the first
step.

## Build the java code:
Activate the conda environment again:
```
conda activate repinpop
```
`RepinPop` requires at least java version 11. Building is done by
[`gradle`](https://gradle.org).

```
$> cd REPIN_ecology/REPIN_ecology
$> gradle build
```

## Set library path.
Some environment variables (in particular `LD_LIBRARY_PATH`) have to be set
explicitely. 

```
source setenv.sh
```

## Launch the server
To launch the server, run

```
$> flask run 
```

And navigate your browser to localhost:5000 .




