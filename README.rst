Data Transfer Kit (DTK)
***********************

The Data Transfer Kit (DTK) is an open-source software library
designed to provide parallel services for solution transfer for
multiphysics simulations. DTK uses a general operator design to
provide scalable parallel services for solution transfer between
shared volumes and surfaces. DTK is supported and used by the
following Department of Energy multiphysics modeling and simulation
programs:

* Consortium for Advanced Simulation of Light Water Reactors (CASL)
  <http://www.casl.gov>

* National Highway Traffic Safety Administration (NHTSA)
  <http://batterysim.org>

* Nuclear Energy Advanced Modeling and Simulation (NEAMS)
  <http://www.ne.anl.gov/NEAMS/>


DataTransferKit Development Team
================================

DTK is developed and maintained by:

* Stuart Slattery <slatterysr@ornl.gov>

* Damien Lebrun-Grandie <lebrungrandt@ornl.gov>

* Roger Pawlowski <rppawlo@sandia.gov>


DataTransferKit Packages
========================

DTK has the following packages:

**Utils**
    General utilities for software development including exception
    handling, MPI-based tools, and functional programming tools

**Interface**
    Core DTK interface package. Interfaces are divided into two
    categories: *Client* and *Operator*. *Client* interfaces define a
    polymorphic API implemented by client applications providing
    access to mesh, geometry, parallel decomposition, shape functions,
    and parametric mappings. *Operator* interfaces define the general
    operator for solution transfer and other concepts for
    constructing solution transfer operators from client code.

**Operators**
    DTK solution transfer operator implementation package. Operators
    contains implementations of the following algorithms:

    * Shape function interpolation
    * Moving least square reconstruction
    * Spline interpolation

**Adapters**
    Client interface adapters for common mesh databases,
    discretization libraries, and geometric objects. Implementations
    include:

    * Sierra Toolkit Mesh (STK Mesh) <http://trilinos.org/packages/stk/>

    * Intrepid: Interoperable Tools for Rapid dEveloPment of
      compatIble Discretizations
      <http://trilinos.org/packages/intrepid/>

    * MOAB: A Mesh-Oriented datABase
      <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>

    * Basic geometric objects

    Outside of the DTK source code, other applications have
    implementations of the DTK client interfaces in their code base
    that may be used to leverage the DTK services and interoperate
    with other libraries that have implemented the interfaces. These
    include:

    * AMP: Advanced Multi-Physics
      <https://rsicc.ornl.gov/codes/ccc/ccc7/ccc-793.html>


Questions, Bug Reporting, and Issue Tracking
============================================

Questions, bug reporting and issue tracking are provided by
GitHub. Please report all bugs by creating a new issue. You can ask
questions by emailing the developers or by creating a issue with the
question tag.


Dependencies
============

DataTransferKit is designed to build and run with a minimum number of
dependencies and is structured largely as a Trilinos package.  The
dependenices and third-party libraries (TPLs) necessary to build DTK
are all open-source and freely available. The Dependencies for DTK are
listed in the following table:

+-----------------------+---------------+---------------------------------------+
| Dependency            | Required      | Comments                              |
+=======================+===============+=======================================+
| C++11                 | Yes           | GNU, Intel, and Clang are suggested   |
+-----------------------+---------------+---------------------------------------+
| DTKData               | Yes           | Large binary files for tests/examples |
+-----------------------+---------------+---------------------------------------+
| TriBITS               | Yes           | Build system provided with Trilinos   |
+-----------------------+---------------+---------------------------------------+
| Trilinos              | Yes           | Current git revision required         |
+-----------------------+---------------+---------------------------------------+
| BLAS/LAPACK           | Yes           | Use vendor-specific implementation    |
+-----------------------+---------------+---------------------------------------+
| MPI                   | No            | OpenMPI and MPICH are suggested       |
+-----------------------+---------------+---------------------------------------+
| MOAB                  | No            | Required to build MOAB adapters       |
+-----------------------+---------------+---------------------------------------+

You can get the most recent Trilinos git revision at
<https://software.sandia.gov/trilinos/repositories/publicTrilinos>. It
is anticipated that when Trilinos 12.0 is released in the spring, DTK
will depend on that release instead of the development
head. Currently, features targeted for the 12.0 Trilinos release that
are not present in the current 11.12 release are required by DTK and
therefore the development head is needed at present to build DTK.

To use build some examples and tests you will need the DTKData
repository which can be found at
<https://github.com/ORNL-CEES/DTKData>. Simply check out the
repository into the top level DataTransferKit directory or provide a
soft link of to the location of the repository.


Building DTK
============

The following steps can be followed to build DTK with MPI support as
well as tests and examples. First, checkout DataTransferKit (assumed
to be in a directory named `DataTransferKit` for these
instructions). Next create a soft link to the cloned copy of the
DTKData repository::

    > cd DataTransferKit
    > ln -s $PATH_TO_DTKDATA

After this, we need to create a soft link of DTK into the main
Trilinos directory. We do this because DTK is a TriBITS package and
will build as a part of the Trilinos build, effectively becoming a
linkable package include among the larger group of Trilinos
packages. We create this link as::

    > cd $PATH_TO_TRILINOS
    > ln -s $PATH_TO_DATATRANSFERKIT

TriBITS is a CMake-based meta-build system
<https://github.com/TriBITSPub> used by Trilinos. Although freely
available on GitHub, a version is also included as a snapshot in
Trilinos and we use that version here. To setup the build we will make
two directories; one for building and one for installing::

    > mkdir $PATH_TO_BUILD_DIR
    > mkdir $PATH_TO_INSTALL_DIR
    > cd $PATH_TO_BUILD_DIR

Next we can run a build shell script that executes CMake with a number
of options to configure both Trilinos and DataTransferKit::

    #!/bin/bash

    # Clear previous configure
    rm -rf CMakeCache.txt
    rm -rf CMakeFiles

    cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=${PATH_TO_INSTALL_DIR} \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D MPI_BASE_DIR:PATH=$PATH_TO_MPI_INSTALL \
    -D TPL_BLAS_LIBRARIES:STRING="${PATH_TO_BLAS_LIB}" \
    -D TPL_LAPACK_LIBRARIES:STRING="${PATH_TO_LAPACK_LIB}" \
    -D Trilinos_ENABLE_CXX11:BOOL=ON \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    -D Trilinos_ASSERT_MISSING_PACKAGES=OFF \
    -D Trilinos_EXTRA_REPOSITORIES="DataTransferKit" \
    -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
    -D DataTransferKit_ENABLE_DBC:BOOL=ON \
    -D DataTransferKit_ENABLE_TESTS:BOOL=ON \
    -D DataTransferKit_ENABLE_EXAMPLES:BOOL=ON \
    $PATH_TO_TRILINOS

Some details on the script: 

* Changing the variable `CMAKE_BUILD_TYPE` from `DEBUG` to `RELEASE`
  will produce an optimized build instead of debug build.

* `MPI_BASE_DIR:PATH=$PATH_TO_MPI_INSTALL` tells CMake where the MPI
  installation you would like to use resides. If you only have one and
  it is set in your environment, simply setting
  `TPL_ENABLE_MPI:BOOL=ON` can be enough.

* `TPL_BLAS_LIBRARIES` and `TPL_LAPACK_LIBRARIES` point to the BLAS and
  LAPACK libraries installed on your system.

* `Trilinos_EXTRA_REPOSITORIES="DataTransferKit"` indicates to TriBITS
  that we are adding DataTransferKit as an additional Trilinos
  package.

* `Trilinos_ENABLE_DataTransferKit:BOOL=ON` instructs TriBITS to build
  DataTransferKit

* `DataTransferKit_ENABLE_DBC` toggles the DataTransferKit
  Design-by-Contract feature `ON` or `OFF`. This feature adds many
  layers of checks into the code useful for debugging purposes at that
  cost of significant additional runtime. Enabling the feature is
  recommended for new users to verify inputs and implementations. This
  feature should be disabled for production calculations once an
  implementation is tested.

* `DataTransferKit_ENABLE_TESTS` toggles if unit tests are `ON` or
  `OFF`

* `DataTransferKit_ENABLE_EXAMPLES=ON` toggles if examples are `ON` or
  `OFF`

To build other packages of DTK, the syntax is similar. For example, to
add the Moab client interface implementations to the build add the
following::

    -D TPL_ENABLE_MOAB:BOOL=ON \
    -D MOAB_LIBRARY_DIRS:PATH=${PATH_TO_MOAB_INSTALL}/lib \
    -D MOAB_INCLUDE_DIRS:PATH=${PATH_TO_MOAB_INSTALL}/include \
    -D Trilinos_ENABLE_DataTransferKitMoabAdapters:BOOL=ON \

Configuration, building, testing, installing (assuming
`DataTransferKit_ENABLE_TESTS` is `ON`) then proceeds as follows
assuming 8 threads are available for building and testing::

    > cd $PATH_TO_BUILD_DIR
    > ./run_cmake_configure.sh
    > make -j8
    > ctest -j8
    > make -j8 install

It is always recommended to build and run unit tests when installing
DTK to ensure that the installation process was correct and that DTK
has no bugs dependent on your system. If your application code using
DTK is not working and unit tests are failing, this can help the
developers track down the problem.
