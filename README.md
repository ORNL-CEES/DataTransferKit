Data Transfer Kit (DTK)
=======================

The Data Transfer Kit (DTK) is a software component designed to
provide parallel services for mesh and geometry searching and data
transfer for arbitrary physics components. In many physics
applications, the concept of mesh and geometry is used to subdivide
the physical domain into a discrete representation to facilitate the
solution of the model problems that describe it. Additionally, the
concept of the field is used to apply degrees of freedom to the mesh
or geometry as a means of function discretization. With the increased
development efforts in multiphysics simulation, adaptive mesh
simulations, and other multiple mesh/geometry problems, generating
parallel topology maps for transferring fields and other data between
meshes is a common operation. DTK is being developed to provide a
suite of concrete algorithm implementations for these services.


DOCUMENTATION
-------------

User and developer documentation is provided by Doxygen
[here](http://cnerg.github.com/DataTransferKit/). A domain model
document can be found
[here](https://github.com/CNERG/DataTransferKit/master/doc/domain_model/domain_model.pdf) for a more formal description.


BUG REPORTING and ISSUE TRACKING
--------------------------------

Bug reporting and issue tracking are provided by GitHub. Please report
all bugs [here](https://github.com/CNERG/DataTransferKit/issues) by
creating a new issue.


DEPENDENCIES
------------

The current Trilinos state is required to build DTK. You can check out
this git public repository
[here](http://trilinos.sandia.gov/publicRepo/index.html)

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
parallel search and transfer operations. See the example directory.


AUTHORS
-------

The following people have made substantial contributions to the
development of DataTransferKit:

* [Stuart Slattery](http://github.com/sslattery)
