#! /bin/bash

cd build
make 
cd ..
source /root/.bashrc
conda activate repinpop
source setenv.sh
cd test/scripts

./dl_zenodo.sh all

neisseria_small=./test.sh neisseria_small clean setup run_java test_java run_andi test_andi run_clustdist test_clustdist plots test_plots
mv /tmp/rarefan_test/neisseria_small/out ../data/neisseria_small/out

dokdonia=./test.sh dokdonia clean setup run_java test_java run_andi test_andi run_clustdist test_clustdist plots test_plots
mv ../data/datasets/dokdonia/out ../data/datasets/dokdonia/ref
mv /tmp/rarefan_test/dokdonia/out ../data/datasets/dokdonia/out

neisseria_plots=./test.sh neisseria clean setup ref_plots test_ref_plots
chlororaphis_plots=./test.sh chlororaphis clean setup ref_plots test_ref_plots

echo "Testing neisseria_small: $neisseria_small"
echo "Testing dokdonia: $dokdonia"
echo "Testing neisseria_plots: $neisseria_plots"
echo "Testing chlororaphis_plots: $chlororaphis_plots"

cd ..
