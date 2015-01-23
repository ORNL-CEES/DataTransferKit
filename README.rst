Data Transfer Kit (DTK)
***********************

The Data Transfer Kit (DTK) is a software library designed to provide
parallel services for solution transfer for arbitrary physics
components. DTK uses a general operator design to provide scalable
services for solution transfer between shared volumes and
services. DTK is both supported and used by the following Department
of Energy multiphysics modeling and simulation programs:

* Consortium for Advanced Simulation of Light Water Reactors (CASL)
  [CASL Homepage](http://www.casl.gov)

* Computer-Aided Engineering of Batteries
  [ORNL CAEBAT Homepage](http://batterysim.org)

* Nuclear Energy Advanced Modeling and Simulation
  [NEAMS Homepage](http://www.ne.anl.gov/NEAMS/)


DataTransferKit Development Team
================================

Profugus is developed and maintained by:

* [Stuart Slattery](https://github.com/sslattery)

* [Damien Lebrun-Grandie](https://github.com/dalg24)

* [Roger Pawlowski](https://github.com/rppawlo)


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
    and parameteric mappings. *Operator* interfaces define the general
    linear operator for solution transfer and other concepts for
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

    * MOAB: A Mesh-Oriented datABase
      [http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB](http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB)

    * Sierra Toolkit Mesh (STK Mesh) [http://trilinos.org/packages/stk/](http://trilinos.org/packages/stk/)

    * Basic geometric objects.

    Outside of the DTK source code, other applications have
    implementations of the DTK client interfaces in their code base
    that may be used to leverage the library services. These include:

    * AMP: Advanced Multi-Physics
      [https://rsicc.ornl.gov/codes/ccc/ccc7/ccc-793.html](https://rsicc.ornl.gov/codes/ccc/ccc7/ccc-793.html)

DOCUMENTATION
-------------

A recent doxygen build is hosted [here](http://ornl-cees.github.io/DataTransferKit/)


BUG REPORTING and ISSUE TRACKING
--------------------------------

Bug reporting and issue tracking are provided by GitHub. Please report
all bugs [here](https://github.com/ORNL-CEES/DataTransferKit/issues) by
creating a new issue.


DEPENDENCIES
------------

To use build some examples and tests you will need the DTKData repository
which can be found [here](https://github.com/ORNL-CEES/DTKData). Simply
check out the repository into the top level DataTransferKit directory.

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

