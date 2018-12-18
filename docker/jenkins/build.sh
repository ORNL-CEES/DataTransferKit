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

if [ "${BUILD_TYPE}" == "clang60-sanitizer" ]
then
    source ../scripts/docker_clang60_env.sh undefined_sanitizer
    ../scripts/docker_cmake -D Trilinos_ENABLE_Fortran=OFF
elif [ "${BUILD_TYPE}" == "gcc71" ]
then
    source ../scripts/docker_gcc71_env.sh
    ../scripts/docker_cmake -D Trilinos_CXX11_FLAGS="-std=c++17" -D Trilinos_ENABLE_OpenMP=OFF -D Tpetra_INST_OPENMP=OFF -D Kokkos_ENABLE_OpenMP=OFF
elif [ "${BUILD_TYPE}" == "gcc54-cuda8" ] || [ "${BUILD_TYPE}" == "gcc54-cuda90" ]
then
    source ../scripts/set_kokkos_env.sh
    ../scripts/docker_cuda_cmake -D Trilinos_ENABLE_SERIAL=OFF -D Tpetra_INST_SERIAL=OFF -D Kokkos_ENABLE_Serial=OFF
elif [ "${BUILD_TYPE}" == "gcc54" ]
then
    ../scripts/docker_cmake -D Trilinos_ENABLE_COVERAGE_TESTING=ON
else
    echo "Unknown BUILD_TYPE"
    exit 1
fi

# build
make -j${NPROC}
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
# upload code coverage only once
if [ "${BUILD_TYPE}" == "gcc54"  ]
then
# collect coverage data
lcov --capture --directory DataTransferKit --output-file lcov.info
# upload it to codecov
curl -s https://codecov.io/bash -o codecov_bash_uploader
chmod +x codecov_bash_uploader
./codecov_bash_uploader -Z -X gcov -f lcov.info
fi
# Run performance tests
if [ "${BUILD_TYPE}" == "gcc54-cuda8" ]
then
  ./DataTransferKit/packages/Search/examples/bvh_driver/DataTransferKit_bvh.exe --benchmark_out=test.json --values=1000000 --queries=100000 --benchmark_repetitions=5 --benchmark_report_aggregates_only=true
  # FIXME Making exit code be zero regardless of the outcome of the performance benchmarking comparison.
  python ../scripts/compare.py -c $GIT_COMMIT -b $BUILD_NUMBER -r ../scripts/ref_gcc54_cuda8.json -n test.json -o ../performance_history.csv || true
fi
