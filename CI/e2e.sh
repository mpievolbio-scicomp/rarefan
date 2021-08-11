#! /bin/bash
# End-to-end testing
# To be run manually inside docker container from top level dir, e.g. docker run mpievolbioscicomp/repinpop_base:latest CI/e2e.sh

OUT=$1

git clone https://gitlab.gwdg.de/mpievolbio-scicomp/repinpop -b testing 2>&1 >> $OUT
cd repinpop
CI/build.sh 2>&1 >> $OUT
CI/test.sh 2>&1 >> $OUT
