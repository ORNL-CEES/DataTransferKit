#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=8}
# cleanup workspace
cd ${TRILINOS_DIR}/DataTransferKit
[ -d build ] && rm -rf build
mkdir build && cd build
# configure trilinos with dtk
../scripts/docker_cmake -D Trilinos_ENABLE_COVERAGE_TESTING=ON
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
# collect coverage data
lcov --capture --directory DataTransferKit --output-file lcov.info
# upload it to codecov
curl -s https://codecov.io/bash -o codecov_bash_uploader
chmod +x codecov_bash_uploader
./codecov_bash_uploader -X gcov -f lcov.info
