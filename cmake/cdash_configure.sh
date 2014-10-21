#!/bin/bash
##---------------------------------------------------------------------------##
## CONFIGURE DTK
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
REQUIRED_VARIABLES=(
TRILINOS_INSTALL_DIR
DTK_BUILD_TYPE
MPI_INSTALL_DIR
BOOST_INSTALL_DIR
DTK_USE_DBC
TRILINOS_SOURCE_DIR
)

for VAR in "${REQUIRED_VARIABLES[@]}"; do
    if [ -z "${!VAR}" ]; then
        echo "Undefined ${VAR} environment variable"
        exit 1
    fi
done

##---------------------------------------------------------------------------##

cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=${TRILINOS_INSTALL_DIR} \
    -D CMAKE_BUILD_TYPE:STRING=${DTK_BUILD_TYPE} \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
    -D CMAKE_CXX_FLAGS:STRING="-Wno-unused-local-typedefs -Wall" \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D Boost_LIBRARY_DIRS:PATH=${BOOST_INSTALL_DIR}/lib \
    -D Boost_INCLUDE_DIRS:PATH=${BOOST_INSTALL_DIR}/include \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    -D Trilinos_EXTRA_REPOSITORIES="DataTransferKit" \
    -D Trilinos_ASSERT_MISSING_PACKAGES=OFF \
    -D Trilinos_ENABLE_CXX11:BOOL=ON \
    -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
    -D DataTransferKit_ENABLE_DBC:BOOL=${DTK_USE_DBC} \
    -D DataTransferKit_ENABLE_TESTS:BOOL=ON \
    -D DataTransferKit_ENABLE_EXAMPLES:BOOL=OFF \
    ${TRILINOS_SOURCE_DIR}

    # -D CMAKE_C_COMPILER=:PATH=${MPI_INSTALL_DIR}/bin/mpicc \
    # -D CMAKE_CXX_COMPILER=:PATH=${MPI_INSTALL_DIR}/bin/mpicxx \
    # -D CMAKE_Fortran_COMPILER=:PATH=${MPI_INSTALL_DIR}/bin/mpif90 \
