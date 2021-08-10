#! /bin/bash
# End-to-end testing
# To be run manually inside docker container from top level dir, e.g. docker run mpievolbioscicomp/repinpop_base:latest CI/e2e.sh

git clone https://gitlab.gwdg.de/mpievolbio-scicomp/repinpop -b testing
cd repinpop
CI/build.sh
CI/test.sh
