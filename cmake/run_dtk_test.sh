#!/bin/bash

export TRILINOS_SOURCE_DIR=/Users/uy7/software/Trilinos
export DTK_SOURCE_DIR=${TRILINOS_SOURCE_DIR}/DataTransferKit
export TRILINOS_INSTALL_DIR=/Users/uy7/builds/DataTransferKit/cdash_dbg
export TRILINOS_BUILD_DIR=${TRILINOS_INSTALL_DIR}
export BOOST_INSTALL_DIR=/Users/uy7/builds/boost_1_55_0
export DTK_BUILD_TYPE=DEBUG
export DTK_USE_DBC=ON
export MPI_INSTALL_DIR=/Users/uy7/builds/mpich-3.1

rm -rf $TRILINOS_BUILD_DIR
mkdir $TRILINOS_BUILD_DIR
cd $TRILINOS_BUILD_DIR
ctest -S $DTK_SOURCE_DIR/cmake/cdash_config.cmake -VV
