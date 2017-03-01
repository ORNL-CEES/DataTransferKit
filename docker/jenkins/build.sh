#!/usr/bin/env bash

set -e

# number of processes with default value
: ${NPROC:=8}
# make a symbolic link to the DTK source dir
ln -s $DTK_DIR ${TRILINOS_DIR}/DataTransferKit
# cleanup workspace
cd ${TRILINOS_DIR}/DataTransferKit
[ -d build ] && rm -rf build
mkdir build && cd build
# configure trilinos with dtk
if [ "${SANITIZE}" == "undefined" ]
then . ../scripts/docker_clang_env.sh undefined_sanitizer
../scripts/docker_cmake -D Trilinos_ENABLE_Fortran=OFF -D TPL_ENABLE_MOAB=OFF
elif [ "${SANITIZE}" == "thread" ]
then . ../scripts/docker_clang_env.sh thread_sanitizer
../scripts/docker_cmake -D Trilinos_ENABLE_Fortran=OFF
else
../scripts/docker_cmake -D Trilinos_ENABLE_COVERAGE_TESTING=ON
fi
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
if [ -z ${SANITIZE} ]
then
# collect coverage data
lcov --capture --directory DataTransferKit --output-file lcov.info
# upload it to codecov
curl -s https://codecov.io/bash -o codecov_bash_uploader
chmod +x codecov_bash_uploader
./codecov_bash_uploader -Z -X gcov -f lcov.info
fi
