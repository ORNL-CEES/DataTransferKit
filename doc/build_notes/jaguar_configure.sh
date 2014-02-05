#!/bin/bash
##---------------------------------------------------------------------------##
## CONFIGURE DataTransferKit on Jaguar
##---------------------------------------------------------------------------##

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

##---------------------------------------------------------------------------##

cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=$PWD \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D MPI_C_COMPILER:FILEPATH="/opt/cray/xt-asyncpe/5.04/bin/cc" \
    -D MPI_CXX_COMPILER:FILEPATH="/opt/cray/xt-asyncpe/5.04/bin/CC" \
    -D MPI_Fortran_COMPILER:FILEPATH="/opt/cray/xt-asyncpe/5.04/bin/ftn" \
    -D BLAS_LIBRARY_DIRS:PATH=/opt/xt-libsci/11.0.04.4/gnu/46/interlagos/lib \
    -D BLAS_LIBRARY_NAMES:STRING="sci_gnu" \
    -D LAPACK_LIBRARY_DIRS:PATH=/opt/xt-libsci/11.0.04.4/gnu/46/interlagos/lib \
    -D LAPACK_LIBRARY_NAMES:STRING="sci_gnu" \
    -D Zoltan_ENABLE_ULLONG_IDS:Bool=ON \
    -D Trilinos_EXTRA_REPOSITORIES="DataTransferKit" \
    -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
    -D DataTransferKit_ENABLE_DBC:BOOL=OFF \
    -D DataTransferKit_ENABLE_TESTS:BOOL=ON \
    -D DataTransferKit_ENABLE_EXAMPLES:BOOL=ON \
    $EXTRA_ARGS \
    /ccs/home/uy7/software/uy7/Trilinos
