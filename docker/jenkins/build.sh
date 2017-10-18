#!/usr/bin/env bash

set -e

# number of processes with default value
: ${NPROC:=8}
# bind mount DTK source dir into Trilinos base dir
mkdir ${TRILINOS_DIR}/DataTransferKit
mount --bind ${DTK_DIR} ${TRILINOS_DIR}/DataTransferKit
# cleanup workspace
cd ${TRILINOS_DIR}/DataTransferKit
[ -d build ] && rm -rf build
mkdir build && cd build
# NOTE: relative paths are invalid after configuration when DTK source dir is
# not directly mounted into Trilinos base source dir.  We build elsewhere and
# move the build directory afterwards...
# configure trilinos with dtk
if [ "${COMPILER}" == "clang" ]
then
    if [ "${SANITIZE}" == "undefined" ]
    then
        source ../scripts/docker_clang_env.sh undefined_sanitizer
        ../scripts/docker_cmake -D Trilinos_ENABLE_Fortran=OFF
    elif [ "${SANITIZE}" == "thread" ]
    then
        source ../scripts/docker_clang_env.sh thread_sanitizer
        ../scripts/docker_cmake -D Trilinos_ENABLE_Fortran=OFF
    else
        source ../scripts/docker_clang_env.sh
        ../scripts/docker_cmake -D Trilinos_ENABLE_Fortran=OFF
    fi
elif [ "${COMPILER}" == "gcc71" ]
then
    source ../scripts/docker_gcc71_env.sh
    ../scripts/docker_cmake -D Trilinos_CXX11_FLAGS="-std=c++17" -D Trilinos_ENABLE_OpenMP=OFF -D Tpetra_INST_OPENMP=OFF -D Kokkos_ENABLE_OpenMP=OFF
else
    if [ -z "${CUDA_VERSION}" ]
    then
        ../scripts/docker_cmake -D Trilinos_ENABLE_COVERAGE_TESTING=ON
    else
        source ../scripts/set_kokkos_env.sh
        ../scripts/docker_cuda_cmake
    fi
fi
# build
make -j${NPROC}
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
# upload code coverage only once
if [ "${COMPILER}" == "gcc54"  ] && [ -z "${CUDA_VERSION}" ]
then
# collect coverage data
lcov --capture --directory DataTransferKit --output-file lcov.info
# upload it to codecov
curl -s https://codecov.io/bash -o codecov_bash_uploader
chmod +x codecov_bash_uploader
./codecov_bash_uploader -Z -X gcov -f lcov.info
fi
