#!/bin/bash
##---------------------------------------------------------------------------##
## CONFIGURE COUPLER ON BEAKER WITH CMAKE
##---------------------------------------------------------------------------##

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

##---------------------------------------------------------------------------##

cmake \
-D CMAKE_INSTALL_PREFIX:PATH=/home/stuart/software/builds/coupler \
-D CMAKE_BUILD_TYPE:STRING=DEBUG \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_Teuchos:BOOL=ON \
-D Teuchos_ENABLE_EXTENDED:BOOL=ON \
-D Trilinos_ENABLE_EpetraExt:BOOL=ON \
-D Trilinos_ENABLE_Tpetra:BOOL=ON \
-D Trilinos_ENABLE_Amesos:BOOL=ON \
-D Trilinos_ENABLE_AztecOO:BOOL=ON \
-D Trilinos_ENABLE_Ifpack:BOOL=ON \
-D Trilinos_ENABLE_ML:BOOL=ON \
-D Trilinos_ENABLE_NOX:BOOL=ON \
-D Trilinos_ENABLE_Anasazi:BOOL=ON \
-D Trilinos_ENABLE_Kokkos:BOOL=ON \
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH=/home/stuart/software/builds/openmpi-1.4.4 \
-D BLAS_LIBRARY_DIRS:PATH=/home/stuart/software/lapack-3.4.0 \
-D BLAS_LIBRARY_NAMES:STRING="blas" \
-D LAPACK_LIBRARY_DIRS:PATH=/home/stuart/software/lapack-3.4.0 \
-D LAPACK_LIBRARY_NAMES:STRING="lapack" \
-D Trilinos_EXTRA_REPOSITORIES="Nemesis;Coupler" \
-D Trilinos_ENABLE_nemesis:BOOL=ON \
-D Trilinos_ENABLE_nemesisHarness:BOOL=ON \
-D Trilinos_ENABLE_nemesisComm:BOOL=ON \
-D Trilinos_ENABLE_coupler:BOOL=ON \
-D coupler_ENABLE_TESTS:BOOL=ON \
$EXTRA_ARGS \
/home/stuart/software/Trilinos