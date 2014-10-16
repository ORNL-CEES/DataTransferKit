Data Transfer Kit (DTK)
=======================

The Data Transfer Kit (DTK) is a software component designed to
provide parallel services for mesh and geometry searching and data
transfer for arbitrary physics components. With the increased
development efforts in multiphysics simulation, adaptive mesh
simulations, and other multiple mesh/geometry problems, generating
parallel topology maps for transferring fields and other data between
meshes and other geometries is a common operation. DTK is being
developed to provide a suite of algorithm implementations for these
services.


DOCUMENTATION
-------------



BUG REPORTING and ISSUE TRACKING
--------------------------------

Bug reporting and issue tracking are provided by GitHub. Please report
all bugs [here](https://github.com/ORNL-CEES/DataTransferKit/issues) by
creating a new issue.


DEPENDENCIES
------------

The current Trilinos state is required to build DTK. You can check out
this git public repository
[here](http://trilinos.sandia.gov/publicRepo/index.html)

C++11 support as well as Boost are required.

DTK can be configured for both serial and parallel builds. For
parallel builds, an MPI implementation is also required. Both OpenMPI
and MPICH have been tested.


CONFIGURE, BUILD and TEST
-------------------------

The DataTransferKit uses the
[TriBITS](http://www.ornl.gov/~8vt/TribitsLifecycleModel_v1.0.pdf)
build system distributed with Trilinos with a required dependency on
[CMake](http://www.cmake.org/). Sample CMake configure scripts can be
found
[here](https://github.com/CNERG/DataTransferKit/master/doc/build_notes)
for various systems. The [CTest](http://www.cmake.org/Wiki/CTest:FAQ)
harness is used for testing.


EXAMPLES
--------

Several examples are provided for using the DataTransferKit for
parallel search and transfer operations. See the examples directory in
each subpackage.


AUTHORS
-------

The following people have made substantial contributions to the
development of DataTransferKit:

* [Stuart Slattery](https://github.com/sslattery)
* [Damien Lebrun-Grandie](https://github.com/dalg24)
* [Roger Pawlowski](https://github.com/rppawlo)
