# Building DTK on Summitdev
The user guide for Summitdev can be found [here](https://www.olcf.ornl.gov/kb_articles/summitdev-quickstart).

## Software environment
Below are the modules that need to be loaded:
```
module swap xl gcc/5.4.0
module load cmake
module load cuda
module load boost
module load netlib-lapack
```

# Configuring Trilinos with DTK
This directory constains a configuration script
([olcf_summitdev_cmake](olcf_summitdev_cmake)) to build DTK and its test suite.
Don't forget to set `CMAKE_INSTALL_PREFIX` and `TRILINOS_DIR`. You also need
to set `nvcc_wrapper` as the underlying C++ compiler for Spectrum MPI:
```
export OMPI_CXX=$TRILINOS_DIR/packages/kokkos/config/nvcc_wrapper
```
