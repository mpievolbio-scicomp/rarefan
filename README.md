# RepinPop

## Installation

## Compilers and build system
The following packages (linux, debian based distro) are required:
* java (version >= 11)
* gcc (or alternativ C compiler)
* libgsl-dev 
* andi
* build-essential
* phyloml

### Install dependencies on debian based linux distros (debian, *ubuntu, mint, ...)
```
sudo apt install linux-libc-dev util-linux git make gcc build-essential libgsl-dev gsl-bin andi wget zip unzip phyloml
```

### Create the conda environment

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

### Build and install cmake targets.
RAREFAN uses cmake to configure, build, and install most of its non-python dependencies including the RAREFAN java code and`clustDist` [Cluster Distances into Phylogenies](https://github.com/EvolBioInf/clustDist.git) for phylogenetic analysis. The installation target directory is the $CONDA_PREFIX directory.

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

### Set library path.
Some environment variables (in particular `LD_LIBRARY_PATH`) have to be set
explicitely. 

```
source setenv.sh
```

## Running RAREFAN from the commandline
The commandline interface to RAREFAN is implemented in *app/utilities/rarefan*. This script can be used to run RAREFAN on a directory that contains genome sequences and rayt protein fasta files.

The syntax of is
```
$> rarefan [-h] [-o OUTDIR] -r REFERENCE [-c MIN_NMER_OCCURRENCE] [-l NMER_LENGTH] -q QUERY_RAYT
               [-e E_VALUE_CUTOFF] [-R] [-j THREADS] [-t TREEFILE] [-i]
               DIR
```
where the commandline arguments are explained as follows:
```
positional arguments:
  DIR                   Contains the genome DNA sequences and RAYT AA sequences to be analyzed.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Results will be written to OUTDIR. OUTDIR will be created if not existing
                        (default: ./rarefan_out).
  -r REFERENCE, --reference REFERENCE
                        Filename of the reference genome sequence
  -c MIN_NMER_OCCURRENCE, --min_nmer_occurrence MIN_NMER_OCCURRENCE
                        Only Nmers of NMER_LENGTH that occur more frequently than MIN_NMER_OCCURRENCE
                        will be taken into account (default: 55). See RAREFAN manual for details.
  -l NMER_LENGTH, --min_nmer_length NMER_LENGTH
                        Only Nmers of NMER_LENGTH that occur more frequently than MIN_NMER_OCCURRENCE
                        will be taken into account (default: 21). See RAREFAN manual for details.)
  -q QUERY_RAYT, --query_rayt QUERY_RAYT
                        Filename or path of the amino acid sequence file containing the RAYT protein
                        sequence (default: None).
  -e E_VALUE_CUTOFF, --e_value_cutoff E_VALUE_CUTOFF
                        e-value cutoff for tblastn of the query rayt sequence against the submitted
                        genomes (default: 1e-30).
  -R, --no-repins       Do not analyse REPINS (default: False).
  -j THREADS, --num_threads THREADS
                        Number of threads for parallel cluster analysis with MCL (default: 24).
  -t TREEFILE, --treefile TREEFILE
                        Filename or path of the phylogenetic tree of submitted genomes (newik format,
                        '.nwk' extension). If none given and more than four genomes are submitted, the
                        tree will be calculated and written to OUTDIR/tmptree.nwk (default:
                        tmptree.nwk).
  -i, --interactive     Interactive mode. Ask for confirmation before starting the analysis run.
```

## Runing the RAREFAN web server
### Database backend
The webserver uses MongoDB as a backend. Install mongodb-server, create a database user named 'rarefan', secured by password, and a database 'rarefan'. Assign the 'dbAdmin' role for the database 'rarefan' to the 'rarefan' user. Consult the [mongodb manual](https://docs.mongodb.com/manual/tutorial/manage-users-and-roles/) if unsure how to do this.

### Configuration
Copy the configuration template *app/config_template.py* to  *app/config.py* and edit the settings. An example is given below.

Jobs submitted to RAREFAN are processed by redis. In your conda environmont, install `rq` and `redis`.

```shell
$> conda install rq redis
```

```python
import os

class Config(object):
    SECRET_KEY = 'supersecretkey'
    SERVER_NAME = 'localhost:5000'
    MONGODB_SETTINGS = {
        'db': 'rarefan',
        'host': 'localhost',
        'port': 27017,
        'username': 'rarefan',
        'password': 'RaReF@npw01'
    }
    REDIS_URL = os.environ.get("REDIS_URL") or 'redis://'
    MAIL_SERVER = 'zimbra.evolbio.mpg.de'
    MAIL_USERNAME='rarefan@evolbio.mpg.de'
    MAIL_PASSWORD='7SaaZv34Xw5isyu'

    MAIL_USE_TLS=True
    MAIL_USE_SSL=False
    MAIL_PORT=25

    MAIL_DEBUG=False
    DEFAULT_MAIL_SENDER='rarefan@evolbio.mpg.de'
```

To launch the server, run

```
$> flask run 
```

And navigate your browser to http://localhost:5000 .

#### NOTE
Data visualisation on a local deployment server is currently not working.

##  Testing
The directory *test/scripts/* contains two scripts:
### *dl_zenodo.sh*
*dl_zenodo.sh* may be used to download reference datasets from zenodo, and unpack the data into the directory *test/data/datasets/*. Datasets can be downloaded individually or together. 

Syntax: 
```shell
./dl_zenodo.sh [all | neisseria | chlororaphis | dokdonia] 
```
### *test.sh*
Having downloaded the reference datasets (see above), the functionality of RAREFAN can be tested using the script *test.sh*. 

Syntax:  
```shell
./test.sh DATASET OPERATION1 [OPERATION2 [OPERATION3 [...] ] ]
```
`DATASET` is the name of one of the reference datasets (*neisseria*, *chlororaphis*, or *dokdonia*). Further test datasets are contained in *test/data* and may also be used, e.g. *neisseria_small* which only contains a small subset of the reference *neisseria* dataset.

`OPERATION` is at least one of:
* `clean`: Remove all data from a previous run
* `setup`: Copy reference data to */tmp/rarefan_test/*`DATASET` and configure the environment for running RAREFAN.
* `run_java`: Run the java code `REPINecology` to compute RAYT and REPIN populations from the input data.
* `run_andi`: Calculate distances between the given genomes.
* `run_clustdist`: Generate a phylogenetic tree based on output from `andi`.
* `plots`: Generate figures from the output data (Phylogenetic trees for input genomes and RAYTs, REPIN and RAYT population sizes, and correlation between replication rate and REPIN population size.)
* `ref_plots`: Generate figures from the reference output data that comes with the downloaded dataset (originally in `DATASET`/*out/*, copied to */tmp/rarefan_test/*`DATASET`*/ref* in the `setup` operation).
* `test_java`: Checks if a file *results.txt* was produced.
* `test_andi`: Checks if the file *.dist* file was generated by `andi`.
* `test_clustdist`: Checks if the *.nwk* treefile was generated by `clustDist`.
* `test_plots`: Checks if all plots were generated.
* `test_ref_plots`: Checks if all plots were generated from the reference output data.
* `test_md5`: Computes md5 checksums for all datafiles  in the output directory (except subdirectories) and compares to checksums in the directory *test/md5/*.


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
docker run --mount source=rarefan_testing,target=/mnt/rarefan_test mpievolbioscicomp/repinpop_base /e2e.sh /mnt/rarefan_test/docker_test.log 
```
This command will launch the docker container, build the RAREFAN code and run a test suite. The entire log output will be written to the file */mnt/rarefan_test/docker_test.log* on the container. Note that */mnt/rarefan_test* is setup as a *persistent docker mount point*, i.e. it can be visited after the test run is complete and the container shut down. The easiest is to navigate to
*/var/lib/docker/volumes/rarefan_testing/_data* on the host machine.







