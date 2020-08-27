#! /bin/sh

TEST_DATA_DIR=data/testdata_small
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

