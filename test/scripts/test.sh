#! /bin/bash
set -e

! getopt --test > /dev/null

if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I am sorry, $(getopt --test) failed in this environment.'
    exit 1
fi

# OPTIONS=d:
# LONGOPTIONS=dataset:

if [ $1 == "-h" ]; then
    echo "Usage: test.sh [DATASET OPERATIONS | -h | --help]"
    exit 1
elif [ $1 == "--help" ]; then
    echo "Usage: test.sh [DATASET OPERATIONS | -h | --help]"
    exit 1
elif [ $# -le 1 ]; then
    echo "Wrong number of arguments ($#). Usage: test.sh [DATASET OPERATIONS | -h | --help]"
    exit 1
fi

# Set root directories 
DATASET="$1"
shift
OPERATIONS="$*"

PROJECT_ROOT_DIR=$(realpath ../../)
TMPDIR=/tmp/rarefan_test

# Test directories
TEST_ROOT_DIR=${PROJECT_ROOT_DIR}/test
TEST_DATA_DIR=${TEST_ROOT_DIR}/data
TEST_MD5_DIR=${TEST_ROOT_DIR}/md5

# Run and output directories.
RUN_DATA_DIR=${TMPDIR}/${DATASET}
RUN_OUT_DIR=${RUN_DATA_DIR}/out
RUN_REF_DIR=${RUN_DATA_DIR}/ref

TREENAME=${DATASET}

echo "DATASET=${DATASET}"
echo "OPERATIONS=${OPERATIONS}"


case $DATASET in 
    chlororaphis)
        ISREF=1
        ;;
    neisseria)
        ISREF=1
        ;;
    dokdonia)
        ISREF=1
        ;;
    *)
        ISREF=0
        ;;
esac


clean() {
    # Clean up previous runs
    rm -rvf ${RUN_DATA_DIR}
    mkdir -vp ${RUN_DATA_DIR}
    mkdir -vp ${RUN_OUT_DIR}
    mkdir -vp ${RUN_REF_DIR}
}

rayt_faa=""
ref_strain=""

setup() {
    # Copy data to run dir.
    if [ $ISREF == "1" ]; then
		TESTCASE_DATA_DIR="${TEST_DATA_DIR}/datasets/${DATASET}"
	    rsync -ruvL ${TESTCASE_DATA_DIR}/in/ ${RUN_DATA_DIR}/
        rsync -ruvL ${TESTCASE_DATA_DIR}/out/ ${RUN_REF_DIR}/
    else
		TESTCASE_DATA_DIR="${TEST_DATA_DIR}/${DATASET}"
		rsync -ruvL ${TESTCASE_DATA_DIR}/ ${RUN_DATA_DIR}/
	fi	
	dataset_vars
		
	echo "RAYT AA sequence will be read from $rayt_faa."
	echo "Reference strain set to $ref_strain."

    rsync -uv $rayt_faa ${RUN_DATA_DIR}
} 

dataset_vars() {
# Dataset specific settings.
	case $DATASET in
		"chlororaphis")
			rayt_faa=${TEST_DATA_DIR}/yafM_SBW25.faa
			ref_strain=chlTAMOak81.fas
			;;
		"dokdonia")
			rayt_faa=${TEST_DATA_DIR}/yafM_Ecoli.faa
			ref_strain=dokd-P16.fas
			;;
		"neisseria")
			rayt_faa=${TEST_DATA_DIR}/yafM_Ecoli.faa
			ref_strain=Nmen_2594.fas
			;;
		*)
			rayt_faa=${TEST_DATA_DIR}/yafM_Ecoli.faa
			ref_strain=Nmen_2594.fas
			;;
	esac
}

java_cmd() {
	dataset_vars
    javacmd="java -Xmx10g -jar ${PROJECT_ROOT_DIR}/REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar ${RUN_DATA_DIR} ${RUN_OUT_DIR} $ref_strain 55 21 $rayt_faa ${TREENAME}.nwk 1e-30 true"
    echo ${javacmd}
}

run_java() {
	cmd=$(java_cmd)
	echo $cmd
	eval $cmd
}

test_java() {
    if test -f ${RUN_OUT_DIR}/results.txt; then
        echo "java results 'results.txt' found."
        exit 0
    else
        echo "java results 'results.txt' not found."
        exit 1
    fi
}
# Run andi and check output.

run_andi() {
    andi ${RUN_DATA_DIR}/*.fas > ${RUN_OUT_DIR}/${TREENAME}.dist
}
test_andi() {
    if test -f ${RUN_OUT_DIR}/${TREENAME}.dist; then
        echo "andi output found."
        exit 0
    else
        echo "andi output not found."
        exit 1
    fi
}

run_clustdist() {
    clustDist ${RUN_OUT_DIR}/${TREENAME}.dist > ${RUN_OUT_DIR}/${TREENAME}.nwk
}

test_clustdist() {
    if test -f ${RUN_OUT_DIR}/${TREENAME}.nwk; then
        echo "clustDist output found."
        exit 0
    else
        echo "clustDist output not found."
        exit 1
    fi
}

test_all_files() {
# Checksums.
    cd ${RUN_OUT_DIR}
    md5sum -c ${TEST_MD5_DIR}/${DATASET}.md5
    cd -
}

ref_plots() {
    if [ $ISREF -ne 1 ]; then
        echo 
        echo "'$DATASET' is not a reference dataset."
        echo 
        exit 1
    fi

    cd ${PROJECT_ROOT_DIR}/shinyapps/analysis
    for rayt in 0 1 2 3 4 5
    do
        Rscript run_analysis.R -d ${RUN_REF_DIR} -t ${TREENAME}.nwk -o ${RUN_REF_DIR}/${DATASET}.png -r ${rayt}
    done
    cd -
}

test_ref_plots() {
    if [ $ISREF -ne 1 ]; then
        echo
        echo "'$DATASET' is not a reference dataset."
        echo
        exit 1
    fi

    for rayt in 0 1 2 3 4 5; do
        for slide in 1 2 3; do
            file=${RUN_REF_DIR}/${DATASET}_rayt${rayt}_00${slide}.png
            if test -f $file; then
                echo "$file found"
            else
                echo "$file not found"
                exit 1
            fi
        done
    done
    exit 0
}

plots() {
    cd ${PROJECT_ROOT_DIR}/shinyapps/analysis
    for rayt in 0 1 2 3 4 5
    do
        Rscript run_analysis.R -d ${RUN_OUT_DIR} -t ${TREENAME}.nwk -o ${RUN_OUT_DIR}/${DATASET}.png -r ${rayt}
    done
    cd -
}

test_plots() {
    for rayt in 0 1 2 3 4 5; do
        for slide in 1 2 3; do
            file=${RUN_OUT_DIR}/${DATASET}_rayt${rayt}_00${slide}.png
            if test -f $file; then
                echo "$file found"
            else
                echo "$file not found"
                exit 1
            fi
        done
    done
    exit 0
}

# Run the desired action
for op in ${OPERATIONS}; do
	$op
done

