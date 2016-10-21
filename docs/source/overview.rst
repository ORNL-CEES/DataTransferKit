Getting started with DTK
========================

Overview
--------

`DataTransferKit <https://github.com/ORNL-CEES/DataTransferKit>`_ is an
open-source software library of parallel solution transfer services for
multiphysics simulations. DTK uses a general operator design to provide
scalable algorithms for solution transfer between shared volumes and surfaces.

DTK was originally developed at the University of Wisconsin - Madison as part of
the Computational Nuclear Engineering Group (`CNERG <http://cnerg.github.io>`_)
and is now actively developed at the Oak Ridge National Laboratory as part of
the Computational Engineering and Energy Sciences (`CEES
<http://energy.ornl.gov>`_) group.

DTK is supported and used by the following projects and programs:

* Oak Ridge National Laboratory (ORNL) Laboratory Directed Research and
  Development (`LDRD
  <https://www.ornl.gov/content/laboratory-directed-research-development>`_)

* Consortium for Advanced Simulation of Light Water Reactors (`CASL
  <http://www.casl.gov>`_)

* Nuclear Energy Advanced Modeling and Simulation (`NEAMS
  <http://www.ne.anl.gov/NEAMS/>`_)

* National Highway Traffic Safety Administration (`NHTSA
  <http://batterysim.org>`_)

DataTransferKit Development Team
--------------------------------

DTK is developed and maintained by:

* `Stuart Slattery <slatterysr@ornl.gov>`_

* `Damien Lebrun-Grandie <lebrungrandt@ornl.gov>`_

* `Bruno Turcksin <turcksinbr@ornl.gov>`_

* `Andrey Prokopenko <prokopenkoav@ornl.gov>`_

* `Roger Pawlowski <rppawlo@sandia.gov>`_

* `Alex McCaskey <mccaskeyaj@ornl.gov>`_


DataTransferKit Packages
------------------------

DTK has the following packages:

**Utils**
    General utilities for software development including exception
    handling, MPI-based tools, and functional programming tools.

**Interface**
    Core DTK interface package. Interfaces are divided into two categories:
    ``Client`` and ``Operator``. ``Client`` interfaces define a polymorphic
    API implemented by client applications providing access to mesh, geometry,
    parallel decomposition, shape functions, and parametric
    mappings. ``Operator`` interfaces define the general operator and vector
    objects for solution transfer and other concepts for constructing solution
    transfer operators from client code.

**Operators**
    DTK solution transfer operator implementation package. Operators
    contains parallel implementations of the following algorithms:

    * L2 Projection
    * Shape function interpolation
    * Moving least square reconstruction with degenerate geometry
      detection
    * Spline interpolation
    * Direct node-to-node transfer

**Adapters**
    Client interface adapters for common mesh databases,
    discretization libraries, and geometric objects. Implementations
    include:

    * Sierra Toolkit Mesh (STK Mesh) `<http://trilinos.org/packages/stk/>`_

    * Intrepid: Interoperable Tools for Rapid dEveloPment of
      compatIble Discretizations
      `<http://trilinos.org/packages/intrepid/>`_

    * MOAB: A Mesh-Oriented datABase
      `<http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_

    * libMesh - A C++ Finite Element Library
      `<http://libmesh.github.io/>`_

    * Basic geometric objects

    * Classic DTK. The classic adapters provide DTK services using the (now
      deprecated) version 1.0 API using the new 2.0 implementation. See the
      notes below for more details.

    * C and Fortran Interfaces. Provides basic C and Fortran interfaces to a
      subset of DTK algorithms.

    Outside of the DTK source code, other applications have
    implementations of the DTK client interfaces in their code base
    that may be used to leverage the DTK services and interoperate
    with other libraries that have implemented the interfaces. These
    include:

    * AMP: Advanced Multi-Physics
      `<https://rsicc.ornl.gov/codes/ccc/ccc7/ccc-793.html>`_

    * Albany mulitphysics code `<https://github.com/gahansen/Albany>`_

**Notes on Classic DTK Adapters**
    This can serve as a starting point for migrating from version 1.0 to 2.0
    but note that the classic adapters only have version 1.0 functionality -
    no version 2.0 functionality is provided. To use these adapters, simply
    build them and your code should link. See the unit tests in the classic
    adapters directory for examples of these changes.

Questions, Bug Reporting, and Issue Tracking
--------------------------------------------

Questions, bug reporting and issue tracking are provided by GitHub. Please
report all bugs by creating a new issue. You can ask questions by creating a
new issue with the question tag.


Publications
------------

Publications to date related to DataTransferKit:

* S. Slattery, *"Mesh-Free Data Transfer Algorithms for Partitioned
  Multiphysics Problems: Conservation, Accuracy, and Parallelism"*, Journal of
  Computational Physics, vol. 307, pp. 164-188, 2016.

* S. Slattery, S. Hamilton, T. Evans, *"A Modified Moving Least Square
  Algorithm for Solution Transfer on a Spacer Grid Surface"*, ANS MC2015 -
  Joint International Conference on Mathematics and Computation (M&C),
  Supercomputing in Nuclear Applications (SNA) and the Monte Carlo (MC)
  Method, Nashville, Tennessee · April 19–23, 2015, on CD-ROM, American
  Nuclear Society, LaGrange Park, IL (2015).

* R. Schmidt, K. Belcourt, R. Hooper, R. Pawlowski, K. Clarno, S. Simunovic, S. Slattery, J. Turner, S. Palmtag,
  *"An Approach for Coupled-Code Multiphysics Core Simulations from a Common
  *Input"*, Annals of Nuclear Energy, Volume 84, pp. 140-152, 2014.

* S. Slattery, P.P.H. Wilson, R. Pawlowski, *“The Data Transfer Kit: A
  Geometric Rendezvous-Based Tool for Multiphysics Data Transfer”*,
  International Conference on Mathematics and Computational Methods Applied to
  Nuclear Science & Engineering (M&C 2013), American Nuclear Society, Sun
  Valley, ID, May 5-9, 2013.
