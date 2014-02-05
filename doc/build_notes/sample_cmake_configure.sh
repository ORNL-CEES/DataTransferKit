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
    -D VERA_EXTRA_REPOSITORIES="Trilinos;DataTransferKit" \
    -D Zoltan_ENABLE_ULLONG_IDS:Bool=ON \
    -D VERA_ENABLE_DataTransferKit:BOOL=ON \
    -D DataTransferKit_ENABLE_DBC:BOOL=ON \
    -D DataTransferKit_ENABLE_TESTS:BOOL=ON \
    -D DataTransferKit_ENABLE_EXAMPLES:BOOL=OFF \
    $EXTRA_ARGS \
    /Users/uy7/software/VERA
