#!/bin/bash
##---------------------------------------------------------------------------##
## CONFIGURE VERA
##---------------------------------------------------------------------------##

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

##---------------------------------------------------------------------------##

cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=$PWD \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D BLAS_LIBRARY_DIRS:PATH=/Users/uy7/builds/lapack-3.4.1 \
    -D BLAS_LIBRARY_NAMES:STRING="blas" \
    -D LAPACK_LIBRARY_NAMES:STRING="lapack" \
    -D LAPACK_LIBRARY_DIRS:PATH=/Users/uy7/builds/lapack-3.4.1 \
    -D Netcdf_LIBRARY_DIRS:PATH=/Users/uy7/builds/netcdf-4.1.3/lib \
    -D Netcdf_INCLUDE_DIRS:PATH=/Users/uy7/builds/netcdf-4.1.3/include \
    -D Boost_INCLUDE_DIRS:PATH=/Users/uy7/software/boost_1_49_0 \
    -D VERA_EXTRA_REPOSITORIES="Trilinos;DataTransferKit" \
    -D VERA_ENABLE_ThreadPool:BOOL=ON \
    -D STK_ENABLE_ThreadPool:BOOL=ON \
    -D VERA_ENABLE_DataTransferKit:BOOL=ON \
    -D DataTransferKit_ENABLE_TESTS:BOOL=ON \
    -D DataTransferKit_ENABLE_EXAMPLES:BOOL=OFF \
    $EXTRA_ARGS \
    /Users/uy7/software/VERA
