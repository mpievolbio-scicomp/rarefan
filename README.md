# RepinPop

## Compilers and build system
The following packages (linux, debian based distro) are required:
* gcc (or alternativ C compiler)
* libgsl-dev 
* andi
* build-essential

### Install dependencies on debian based linux distros (debian, *ubuntu, mint, ...)
```
sudo apt install linux-libc-dev util-linux git make gcc build-essential libgsl-dev gsl-bin andi wget zip unzip
```

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
Not all dependencies are available on the conda archives. `clustDist` [Cluster Distances into Phylogenies](https://github.com/EvolBioInf/clustDist.git)
is handled by a the script `CMakeLists.txt` to be consumed by the `cmake` utility.


We observed that `clustDist` produces faulty results if compiled inside the conda environment. As a workaround, we recommend to build `clustDist` with deactivated conda environment. Here, we If you wish to install into Record the value of the `$CONDA_PREFIX` environment variable, e.g.

```
echo $CONDA_PREFIX > conda_prefix.txt
```

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

## Testing
### Docker
We provide a docker container that packs all dependencies of the java backengine (java code).

### Pull the container
To pull the most recent docker container, run (in a terminal)
```
docker pull mpievolbioscicomp/repinpop_base
```

### Build the application in the docker container Run the test
```
docker run -i -t





