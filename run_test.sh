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

TEST_DATA_DIR=data/testdata_${dataset}
RUN_DATA_DIR=/tmp/testdata
RUN_OUT_DIR=$RUN_DATA_DIR/out

cp -rv $TEST_DATA_DIR $RUN_DATA_DIR
cp -v data/yafM_Ecoli.faa $RUN_DATA_DIR

status=1
java -jar REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar\
    $RUN_DATA_DIR\
    $RUN_OUT_DIR\
    Nmen_2594.fas\
    55 21\
    $RUN_DATA_DIR/yafM_Ecoli.faa\
    tmptree.nwk\
    1e-30\
    false &&\
Rscript ./displayREPINsAndRAYTs.R $RUN_OUT_DIR &&\
display $RUN_OUT_DIR/repins.png &&\
display $RUN_OUT_DIR/correlations.png &&\
status=0

# rm -rf $RUN_DATA_DIR
exit ${status}
