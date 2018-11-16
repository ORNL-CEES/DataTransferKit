# Building DTK on CADES SHPC Condo

The CADES SHPC Condo user guide in available [here](http://support.cades.ornl.gov/user-documentation/_book/#).

## Software stack
From the available modules, below are the ones that need to be loaded:
```bash
module load git cmake
module load PE-gnu/2.0
module load cuda/9.2
module load boost
module load ATLAS
```
It is recommended to create a script that you will be able to source later to load these modules.
You probably want to export the environment variable `TRILINOS_DIR` to specify the path to your Trilinos project source directory in that same script.

## Download Trilinos/DTK, configure, build
This directory contains a configuration script ([cades_shpc_condo_cmake](cades_shpc_condo_cmake)) to build on Condo.
Here are some instructions you could run to build DTK:
```bash
git clone https://github.com/trilinos/Trilinos.git
cd Trilinos
export TRILINOS_DIR=$PWD
git clone https://github.com/ORNL-CEES/DataTransferKit.git
cd DataTransferKit
mkdir build && cd build
source ../scripts/set_kokkos_env.sh
../script/cades_shpc_condo_cmake
make -j4
```

## Execute a job on your SHPC Condo allocation to run the tests/examples
Please refer to the [CADES user guide](http://support.cades.ornl.gov/user-documentation/_book/condos/how-to-use/execute-a-job.html)
for more details on how to create a portable batch system (PBS) script and execute a job on SHPC Condo.
Below is how you may start an interactive session on a node that has GPUs:
```bash
qsub -I -l nodes=1:ppn=16 -l walltime=0:10:00 -W group_list=cades-ccsd -A ccsd -q gpu_p100
```
Then you would need to load appropriate programming environement again using `module load` and you would do:
```bash
cd $TRILINOS_DIR/DataTransferKit/build
ctest
```
