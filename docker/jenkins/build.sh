#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=8}
# set a few environment variables
TRILINOS_VERSION=12.8.1
DTK_SOURCE_DIR=${PREFIX}/source/dtk
TRILINOS_BUILD_DIR=${DTK_SOURCE_DIR}/build
export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_VERSION}
export TRILINOS_INSTALL_DIR=${PREFIX}/install/trilinos/${TRILINOS_VERSION}
# cleanup workspace
[[ -d ${TRILINOS_BUILD_DIR} ]] && rm -rf ${TRILINOS_BUILD_DIR}
mkdir ${TRILINOS_BUILD_DIR} && cd ${TRILINOS_BUILD_DIR}
# configure trilinos with dtk
${TRILINOS_SOURCE_DIR}/DataTransferKit/scripts/docker_configure_cmake.sh
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
