#! /bin/bash

cd build
make 
cd ..
source /root/.bashrc
conda activate repinpop
source setenv.sh
cd test/scripts

./dl_zenodo.sh all

# Run all.
./test.sh neisseria_small clean setup run_java run_andi run_clustdist plots
./test.sh neisseria clean setup run_java run_andi run_clustdist plots
./test.sh chlororaphis clean setup run_java run_andi run_clustdist plots
./test.sh dokdonia clean setup run_java run_andi run_clustdist plots

# Test all.
./test.sh neisseria_small test_java test_andi test_clustdist test_plots test_md5
./test.sh neisseria test_java test_andi test_clustdist test_plots test_md5
./test.sh chlororaphis test_java test_andi test_clustdist test_plots test_md5
./test.sh dokdonia test_java test_andi test_clustdist test_plots test_md5

cd ..
