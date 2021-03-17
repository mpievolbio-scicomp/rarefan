#!/bin/bash

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=60:00
#SBATCH --output=outfile-%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c.fortmanngrote
#SBATCH -C scratch

module load conda/2019.03
conda activate repinpop
source ~/Repositories/repinpop/setenv.sh


SCRIPT="java -Xmx14g -jar REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar    gwdu102:/scratch/rarefan/rarefan_test    gwdu102:/scratch/rarefan/rarefan_test/out    Nmen_2594.fas    55 21    gwdu102:/scratch/rarefan/rarefan_test/yafM_Ecoli.faa    tmptree.nwk    1e-30    false &&
Rscript ./displayREPINsAndRAYTs.R gwdu102:/scratch/rarefan/rarefan_test/out &&
