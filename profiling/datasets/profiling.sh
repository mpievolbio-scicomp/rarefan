#! /bin/bash -l


# for ds in chlor smalt neiss; do
for ds in smalt neiss; do
    echo $ds
    
    # for replicate in 1 2 3; do
        for number_of_genomes in 5 10 20 40; do
            echo "Picking $number_of_genomes random genomes from $ds."

            # for nthreads in 4 8 12 24; do
            nthreads=4
                for q in ../../../test/data/yafM_Ecoli.faa ../../../test/data/yafM_SBW25.faa; do

                    mkdir tmp
                    cd tmp
                    for fas in $(cat ../${ds}/genomes.txt | shuf -n $number_of_genomes); do
                        ln -sf $fas
                    done
                    
                    ref=$(ls *.fas | shuf -n 1)


                    # run rarafan
                    export PYTHONPATH=../../../app/utilities:$PYTHONPATH
                    
                    # Start timer
                    t0=$(date -Ins)
                    python -m rarefan_cli -r $ref -q $q -j $nthreads $PWD
                    t1=$(date -Ins)

                    # Go up to prepare for next run.
                    cd .. 
                    
                    # rm run data.
                    rm -r tmp
                        
                    # Log parameters and timestamps.
                    echo "${ds}\t${replicate}\t${number_of_genomes}\t${nthreads}\t${q}\t${ref}\t${t0}\t${t1}" | tee -a profiles.tsv
                    
                done
            # done
        done
    # done
done
