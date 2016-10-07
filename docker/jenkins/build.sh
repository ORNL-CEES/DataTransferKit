#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=8}
# cleanup workspace
cd ${TRILINOS_DIR}/DataTransferKit
[ -d build ] && rm -rf build
mkdir build && cd build
# configure trilinos with dtk
../scripts/docker_cmake
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
