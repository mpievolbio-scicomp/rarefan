#! /bin/bash
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

# Parse command line arguments: First argument is the dataset, subsequent arguments are operations.
DATASET="$1"
shift
OPERATIONS="$*"

# Set root directories
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
    ecoli)
        ISREF=1
        ;;
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
mem="10g"

setup() {
    # Copy data to run dir.
    if [ $ISREF == "1" ]; then
		TESTCASE_DATA_DIR="${TEST_DATA_DIR}/datasets/${DATASET}"
	    rsync -ruL ${TESTCASE_DATA_DIR}/in/ ${RUN_DATA_DIR}/
        rsync -ruL ${TESTCASE_DATA_DIR}/out/ ${RUN_REF_DIR}/
    else
		TESTCASE_DATA_DIR="${TEST_DATA_DIR}/${DATASET}"
		rsync -ruL ${TESTCASE_DATA_DIR}/ ${RUN_DATA_DIR}/
	fi	
	dataset_vars
		
	echo "RAYT AA sequence will be read from $rayt_faa."
	echo "Reference strain set to $ref_strain."

    rsync -u $rayt_faa ${RUN_DATA_DIR}
} 

dataset_vars() {
# Dataset specific settings.
	case $DATASET in
		"ecoli")
			rayt_faa=${TEST_DATA_DIR}/yafM_Ecoli.faa
			ref_strain=MG1655.fas
			analyze_repins='false'
			mem="4g"
			;;
    "chlororaphis")
			rayt_faa=${TEST_DATA_DIR}/yafM_SBW25.faa
			ref_strain=chlTAMOak81.fas
			analyze_repins='true'
			;;
		"dokdonia")
			rayt_faa=${TEST_DATA_DIR}/yafM_Ecoli.faa
			ref_strain=dokd-P16.fas
			mem="4g"
			analyze_repins='true'
			;;
		"neisseria")
			rayt_faa=${TEST_DATA_DIR}/yafM_Ecoli.faa
			ref_strain=Nmen_2594.fas
			analyze_repins='true'
			;;
		*)
			rayt_faa=${TEST_DATA_DIR}/yafM_Ecoli.faa
			ref_strain=Nmen_2594.fas
			mem="4g"
			analyze_repins='true'
			;;
	esac
}

java_cmd() {
	dataset_vars
    javacmd="java -Xmx${mem} -jar ${PROJECT_ROOT_DIR}/REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar \
                                  ${RUN_DATA_DIR}\
                                  ${RUN_OUT_DIR}\
                                  ${ref_strain}\
                                  55 21\
                                  ${rayt_faa}\
                                  ${TREENAME}.nwk 1e-30 ${analyze_repins} 1"
    echo ${javacmd}
}

run_java() {
	cmd=$(java_cmd)
	echo $cmd
	eval $cmd
}

test_java() {
    echo
    echo "*** Checking 'REPIN_ecology' output. ***"
    if test -f ${RUN_OUT_DIR}/results.txt; then
        echo "java results 'results.txt' found."
    else
        echo "java results 'results.txt' not found."
        exit 1
    fi
}
# Run andi and check output.

run_andi() {
    andi -j ${RUN_DATA_DIR}/*.fas ${RUN_DATA_DIR}/*.fn ${RUN_DATA_DIR}/*.fna ${RUN_DATA_DIR}/*.fastn ${RUN_DATA_DIR}/*.fasta> ${RUN_OUT_DIR}/${TREENAME}.dist
}
test_andi() {
    echo
    echo "*** Checking 'andi' output. ***"
    if test -f ${RUN_OUT_DIR}/${TREENAME}.dist; then
        echo "andi output found."
    else
        echo "andi output not found."
        exit 1
    fi
}

run_clustdist() {
    clustDist ${RUN_OUT_DIR}/${TREENAME}.dist > ${RUN_OUT_DIR}/${TREENAME}.nwk
}

test_clustdist() {
    echo
    echo "*** Checking 'clustDist' output. ***"
    if test -f ${RUN_OUT_DIR}/${TREENAME}.nwk; then
        echo "clustDist output found."
    else
        echo "clustDist output not found."
        exit 1
    fi
}

test_md5() {
# Checksums.
    echo
    echo "*** Checking MD5 checksums ***"
    cd ${RUN_OUT_DIR}

    if md5sum -c ${TEST_MD5_DIR}/${DATASET}.md5; then
      echo
      echo "==> All MD5 checks passed <=="
      echo
    else
      echo
      echo "==> MD5 check failed. <=="
      echo
      exit 1
    fi
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
    echo
    echo "*** Checking Reference plots. ***"
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
}

# Run the desired action
for op in ${OPERATIONS}; do
  $op
done