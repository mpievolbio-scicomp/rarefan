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

### Issues when building `clustDist` inside the conda env.
We observed that `clustDist` produces faulty results if compiled inside the conda environment. As a workaround, we recommend to build `clustDist` with deactivated conda environment. Nevertheless, we *install* `clustDist` and other build products into the conda environment. 

Record the value of the `$CONDA_PREFIX` environment variable, e.g.

```
echo $CONDA_PREFIX > conda_prefix.txt
```

### Deactivate the conda environment.
```
conda deactivate
```

### Build in separate build directory.
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
Change back into the project's root directory
```shell
cd ..
```
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

Our gitlab-CI configuration is setup to use this docker image as a basis of the build procedure.

### End-to-end testing in docker

Scripts in `CI/` are used for testing and building. You can also utilize these scripts for local testing:

```
docker run --mount source=tmp,target=/mnt mpievolbioscicomp/repinpop_base /e2e.sh 2>&1 | tee /mnt/rarefan_test/docker_test.log 
```
This command will launch the docker container, build the RAREFAN code and run a test suite. The entire log output will be written to the file */mnt/rarefan_test/docker_test.log* on the container. Note that */mnt/rarefan_test* is setup as a *persistent docker mount point*, i.e. it can be visited after the test run is complete and the container shut down, e.g. in an interactive docker run
```shell
docker run -i -t --mount source=tmp,target=/mnt
```

or on the host machine under */var/lib/docker/volumes/tmp/_data*.






