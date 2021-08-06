#! /bin/bash
set -e

! getopt --test > /dev/null

if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I am sorry, $(getopt --test) failed in this environment.'
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

TEST_DATA_DIR=test/data/${dataset}
RUN_DATA_DIR=/tmp/rarefan_test
RUN_OUT_DIR=${RUN_DATA_DIR}/out

rsync -ruvL ${TEST_DATA_DIR}/ ${RUN_DATA_DIR}/
rsync -uv test/data/yafM_Ecoli.faa ${RUN_DATA_DIR}
rsync -uv test/data/yafM_SBW25.faa ${RUN_DATA_DIR}


if [[ ! $dataset == *"chlororaphis"* ]]; then
	treename="tmptree"
	java -Xmx10g -jar REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar ${RUN_DATA_DIR} ${RUN_OUT_DIR} Nmen_2594.fas 55 21 ${RUN_DATA_DIR}/yafM_Ecoli.faa ${treename}.nwk 1e-30 true
fi
if [[ $dataset == *"chlororaphis"* ]]; then
    treename="chlororaphis"
    RUN_OUT_DIR=${RUN_DATA_DIR}/test_out
    java -Xmx10g -jar REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar ${RUN_DATA_DIR}/in ${RUN_OUT_DIR} chlTAMOak81.fas 55 21 ${RUN_DATA_DIR}/yafM_SBW25.faa ${treename}.nwk 1e-30 true
fi

andi ${RUN_DATA_DIR}/*.fas > ${RUN_OUT_DIR}/${treename}.dist
clustDist ${RUN_OUT_DIR}/${treename}.dist > ${RUN_OUT_DIR}/${treename}.nwk