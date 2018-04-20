# set nvcc as the underlying C++ compiler for Open MPI
export OMPI_CXX=$TRILINOS_DIR/packages/kokkos/bin/nvcc_wrapper

# set nvcc as the underlying C++ compiler for MPICH
export MPICH_CXX=$TRILINOS_DIR/packages/kokkos/bin/nvcc_wrapper

# set the original C++ compiler as the underlying compiler for the nvcc
# wrapper
export NVCC_WRAPPER_DEFAULT_COMPILER=$PATH_TO_CXX_COMPILER
