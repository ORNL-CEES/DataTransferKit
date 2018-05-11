# Building DTK on Summitdev
The user guide for Summitdev can be found [here](https://www.olcf.ornl.gov/kb_articles/summitdev-quickstart).

## Software Environment
You need to load the following modules:
```
module load gcc/5.4.0 cmake cuda boost netlib-lapack
```

## Download, Configure, Build Trilinos+DTK
```bash
git clone https://github.com/trilinos/Trilinos.git
cd Trilinos
export TRILINOS_DIR=$PWD

git clone https://github.com/ORNL-CEES/DataTransferKit.git
cd DataTransferKit

mkdir build
cd build

module load cmake gcc/5.4.0 cuda boost netlib-lapack
export OMPI_CXX=$TRILINOS_DIR/packages/kokkos/bin/nvcc_wrapper

../scripts/olcf_summitdev_cmake
make -j8
```

## Run Tests

In `Trilinos/DataTransferKit/build` run:
```
bsub -P YOUR_PROJECT_ID -nnodes 1 -W 120 -Is $SHELL
module load cmake gcc/5.4.0 cuda boost netlib-lapack
export OMP_NUM_THREADS=5 # jsrun needs this

ctest
```
