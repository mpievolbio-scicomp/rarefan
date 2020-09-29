#! /bin/bash
set -o errexit -o pipefail -o noclobber -o nounset

! getopt --test > /dev/null

if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I’m sorry, `getopt --test` failed in this environment.'
    exit 1
fi

OPTIONS=d:
LONGOPTIONS=dataset:

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi

d=-
if [ $# -lt 2 ]; then
    echo "Usage: run_test.sh [-d | --dataset] DATASET [-h | --help]"
    exit 1
fi
while [[ "$#" -gt 0 ]]; do
    case $1 in 
        -d|--dataset)
            dataset="$2"
            shift
            ;;
        -h|--help)
            echo "Usage: run_test.sh [-d | --dataset] DATASET [-h | --help]"
            exit 1
            ;;
        *)
            echo "Unknown parameter passed to $1"
            exit 3
            ;;
    esac
    shift
done

HPC_SUBMIT_NODE=gwdu101
TEST_DATA_DIR=data/neisseria_${dataset}
RUN_DATA_DIR=${HPC_SUBMIT_NODE}:/scratch/rarefan/rarefan_test
RUN_OUT_DIR=$RUN_DATA_DIR/out

rsync -ruvL  $TEST_DATA_DIR/ $RUN_DATA_DIR/
rsync -uv data/yafM_Ecoli.faa $RUN_DATA_DIR

cat > script.sh << EOF
#!/bin/bash

#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=60:00
#SBATCH --output=outfile-%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c.fortmanngrote
#SBATCH -C scratch

module load conda/2019.03
source activate repinpop
source ~/Repositories/repinpop/setenv.sh


SCRIPT="java -Xmx14g -jar REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar\
    $RUN_DATA_DIR\
    $RUN_OUT_DIR\
    Nmen_2594.fas\
    55 21\
    $RUN_DATA_DIR/yafM_Ecoli.faa\
    tmptree.nwk\
    1e-30\
    false &&\

Rscript ./displayREPINsAndRAYTs.R $RUN_OUT_DIR &&\

EOF

scp script.sh $RUN_DATA_DIR

#ssh ${HPC_SUBMIT_NODE} sbatch /scratch/rarefan/rarefan_test/script.sh
