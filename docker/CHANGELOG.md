# Changelog

# 20.03.0 (March 23rd, 2020)
### Change
- Updated image to nvidia/cuda:10.1-devel
- Remove vim and emacs
- Updated CMake to 3.16.4
- Updated OpenMPI to 4.0.2
- Updated HDF5 to 1.10.6
- Updated NetCDF to 4.7.3
- Added Kokkos 3.0
- Updated Trilinos to 2a24058e to 7d5bf250 (merge of Kokkos 3.0.0)
- Updated google benchmark to 1.5.0

## 19.06.0 (June 4th, 2019)
### Change
- Updated Kokkos Profiling to 2.7.24
- Updated Trilinos from a27130b to 2a24058e (merge of Kokkos 2.8.0)
- Updated netcdf repo location
  NOTE: netcdf's ftp server is unaccessible

## 19.01.0 (January 14th, 2019)
### Change
- Upgrade CMake to v3.13.3
- Upgrade all LLVM packages to 6.0.1
- Upgrade Open MPI to 2.1.6
NOTE: made base image a build argument with default value ubuntu:16.04
```bash
# build the new images and push them to Docker Hub
docker build --pull -t dalg24/dtk-stack:19.01.0 -f Dockerfile_stack .
docker build --pull --build-arg BASE=nvidia/cuda:8.0-devel-ubuntu16.04 -t dalg24/dtk-stack:19.01.0-cuda80 -f Dockerfile_stack .
docker build --pull --build-arg BASE=nvidia/cuda:9.0-devel-ubuntu16.04 -t dalg24/dtk-stack:19.01.0-cuda90 -f Dockerfile_stack .
docker push dalg24/dtk-stack:19.01.0
docker push dalg24/dtk-stack:19.01.0-cuda80
docker push dalg24/dtk-stack:19.01.0-cuda90
# after testing and approval tag them as latest and push to Docker Hub
docker tag dalg24/dtk-stack:19.01.0 dalg24/dtk-stack:latest
docker tag dalg24/dtk-stack:19.01.0-cuda80 dalg24/dtk-stack:latest-cuda80
docker tag dalg24/dtk-stack:19.01.0-cuda90 dalg24/dtk-stack:latest-cuda90
docker push dalg24/dtk-stack:latest
docker push dalg24/dtk-stack:latest-cuda80
docker push dalg24/dtk-stack:latest-cuda90
```

## 18.09.0 (September 5th, 2018)
### TPLs addition
- Add benchmark support library from google (version 1.4.1)

## 18.04.0 (April 18th, 2018)
### TPLs update
- Updated Boost from 1.65.1 to 1.67.0
- Updated Trilinos from 89b8c7f (12.12.1) to a27130b (HEAD on master)
- Updated Kokkos Profiling and Debugging Tools to version tagged as release 2.5.00 on GitHub
- Updated HDF5 from 1.10.1 to 1.10.2
- Updated NetCDF from 4.4.1.1 to 4.6.1
### Change
- Upgrade CMake to v3.11.1
- Upgrade all LLVM packages to 6.0.0
COMMENT: downloaded binaries for Clang + LLVM (including support for OpenMP) instead of using system package manager/building from source
- Upgrade Open MPI to 2.1.3

## 17.11.0 (November 16th, 2017)
### Change
- Fix Boost installation to support RPATH
- Enable auto formatting for Vim within the dev container

## 17.10.1 (October 16th, 2017)
### Change
- Add Open MPI 2.1.2 to replace version 1.10.1 from the system package manager
COMMENT: intent is to be able to take advantage of MPI ABI compatibility on HPC ressources
### TPLs removal
- Remove PETSc and libMesh

## 17.10.0 (October 5th, 2017)
### TPLs update
- Updated Boost from 1.63.0 to 1.65.1
- Updated Trilinos from 4c8adc9 to 89b8c7f which has been tagged as release 12.12.1 on GitHub
- Updated Kokkos Profiling and Debugging Tools to latest version
- Updated PETSc from 3.6.4 to 3.7.5
- Updated libMesh from 1.2.0 to 1.2.1
COMMENT: decided to keep PETSc/libMesh for now and updated them to their recommended versions in MOOSE installation guidelines
### Change
- Upgrade versions of LLVM, Clang, and OpenMP runtime to 5.0.0
- Use ClangFormat version 5.0 to enforce code formatting style (instead of version 4.0)
