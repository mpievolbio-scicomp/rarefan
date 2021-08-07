#! /bin/bash
set -e

# Set root directories 
DATASET="chlororaphis"
PROJECT_ROOT_DIR=../..
TMPDIR=/tmp/rarefan_test

# Test directories
TEST_ROOT_DIR=${PROJECT_ROOT_DIR}/test
TEST_DATA_DIR=${TEST_ROOT_DIR}/data
TEST_MD5_DIR=${TEST_ROOT_DIR}/md5
TESTCASE_DATA_DIR=${TEST_DATA_DIR}/datasets/${DATASET}/in

# Run and output directories.
RUN_DATA_DIR=${TMPDIR}/${DATASET}
RUN_OUT_DIR=${RUN_DATA_DIR}/out

# Clean up previous runs
rm -rvf ${RUN_DATA_DIR}
mkdir -vp ${RUN_DATA_DIR}

# Copy data to run dir.
rsync -ruvL ${TESTCASE_DATA_DIR}/ ${RUN_DATA_DIR}/
rsync -uv ${TEST_DATA_DIR}/yafM_SBW25.faa ${RUN_DATA_DIR}

TREENAME=${DATASET}

javacmd="java -Xmx10g -jar ${PROJECT_ROOT_DIR}/REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar ${RUN_DATA_DIR} ${RUN_OUT_DIR} chlTAMOak81.fas 55 21 ${RUN_DATA_DIR}/yafM_SBW25.faa ${TREENAME}.nwk 1e-30 true"
echo ${javacmd}

java -Xmx10g -jar ${PROJECT_ROOT_DIR}/REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar ${RUN_DATA_DIR} ${RUN_OUT_DIR} chlTAMOak81.fas 55 21 ${RUN_DATA_DIR}/yafM_SBW25.faa ${TREENAME}.nwk 1e-30 true

# Run andi and check output.
andi ${RUN_DATA_DIR}/*.fas > ${RUN_OUT_DIR}/${TREENAME}.dist
test -f ${RUN_DATA_DIR}/${TREENAME}.dist

clustDist ${RUN_OUT_DIR}/${TREENAME}.dist > ${RUN_OUT_DIR}/${TREENAME}.nwk
test -f ${RUN_DATA_DIR}/${TREENAME}.nwk

# Checksums.
cd ${RUN_OUT_DIR}
md5sum -c ${TEST_MD5_DIR}/${DATASET}.md5
cd -

# Make plots
cd ${RUN_OUT_DIR}

for rayt in 0 1 2 3 4 5
do
	Rscript ${PROJECT_ROOT_DIR}/shinyapps/analysis/run_analysis.R -d ${RUN_OUT_DIR} -t ${TREENAME}.nwk -o ${DATASET}.png -r ${rayt}
done

test -f ${DATASET}_0_001.png
test -f ${DATASET}_0_002.png
test -f ${DATASET}_0_003.png
test -f ${DATASET}_1_001.png
test -f ${DATASET}_1_002.png
test -f ${DATASET}_1_003.png
test -f ${DATASET}_2_001.png
test -f ${DATASET}_2_002.png
test -f ${DATASET}_2_003.png
test -f ${DATASET}_3_001.png
test -f ${DATASET}_3_002.png
test -f ${DATASET}_3_003.png
test -f ${DATASET}_4_001.png
test -f ${DATASET}_4_002.png
test -f ${DATASET}_4_003.png
test -f ${DATASET}_5_001.png
test -f ${DATASET}_5_002.png
test -f ${DATASET}_5_003.png

cd -
