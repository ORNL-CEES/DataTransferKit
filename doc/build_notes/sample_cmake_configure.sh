#!/bin/bash
##---------------------------------------------------------------------------##
## CONFIGURE COUPLER ON BEAKER WITH CMAKE
##---------------------------------------------------------------------------##

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

##---------------------------------------------------------------------------##

cmake \
-D CMAKE_INSTALL_PREFIX:PATH=/home/stuart/software/builds/VERA \
-D CMAKE_BUILD_TYPE:STRING=DEBUG \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH=/home/stuart/software/builds/openmpi-1.4.4 \
-D BLAS_LIBRARY_DIRS:PATH=/home/stuart/software/lapack-3.4.0 \
-D BLAS_LIBRARY_NAMES:STRING="blas" \
-D LAPACK_LIBRARY_DIRS:PATH=/home/stuart/software/lapack-3.4.0 \
-D LAPACK_LIBRARY_NAMES:STRING="lapack" \
-D VERA_EXTRA_REPOSITORIES="Trilinos;Coupler" \
-D Teuchos_ENABLE_EXTENDED:BOOL=ON \
-D VERA_ENABLE_Coupler:BOOL=ON \
-D Coupler_ENABLE_TESTS:BOOL=ON \
-D Coupler_ENABLE_EXAMPLES:BOOL=ON \
$EXTRA_ARGS \
/home/stuart/software/VERA
