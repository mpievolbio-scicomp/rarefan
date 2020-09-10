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

## Install 3rd party libraries through cmake.
Not all dependencies are available on the conda archives. `andi` [Efficient
Estimation of Evolutionary Distances](https://github.com/EvolBioInf/andi.git),
it's dependency `divsufsort` [A lightweight suffix-sorting library](https://github.com/y-256/libdivsufsort.git), and `clustDist` [Cluster Distances into Phylogenies](https://github.com/EvolBioInf/clustDist.git) are handled by a the script `CMakeLists.txt` to be consumed by the `cmake` utility:

```
$> mkdir build
$> cd build
$> cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
$> make
```

This will download the required source codes for all three dependencies, build,
and install the executables into the `conda` environment created in the first
step.

## Build the java code:
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
Our application is served as a web form where users can upload their sequence
files, set the parameters, run the analysis and visualize the results. The web
form is implemented in a jupyter notebook and served with the `voila` framework.
To launch the server, run

```
$> voila --enable_nbextensions --port=<YOUR_FREE_PORT> --no-browser --Voila.ip=0.0.0.0
```

NOTE: You have to specify the port on which to run the service. The last
argument will make the server listen to requests coming from any IP. Restrict
the range of IPs or leave out this argument if this behaviour is unwanted. For
more advanced server configurations see the [`voila`
documentation](https://voila.readthedocs.io/en/stable/deploy.html).




