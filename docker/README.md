# Changelog

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

## 17.10.1 (October 16th, 2017)
### Change
- Add Open MPI 2.1.2 to replace version 1.10.1 from the system package manager
COMMENT: intent is to be able to take advantage of MPI ABI compatibility on HPC ressources
### TPLs removal
- Remove PETSc and libMesh
## TODO
Clear things up about install prefix of LLVM OpenMP

## 17.11.0 (November 16th, 2017)
### Change
- Fix Boost installation to support RPATH
- Enable auto formatting for Vim within the dev container
