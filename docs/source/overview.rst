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

* Exascale Computing Project (`ECP ALExa
  <https://www.exascaleproject.org/project/alexa-accelerated-libraries-exascale>`_)

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

Alumni:

* `Roger Pawlowski <rppawlo@sandia.gov>`_

* `Alex McCaskey <mccaskeyaj@ornl.gov>`_


DataTransferKit Packages
------------------------

DTK has the following packages:

**Utils**
    General utilities for software development including exception handling,
    and other functional programming tools

**Meshfree**
    Point cloud based operators (e.g., nearest neighbor, moving least squares,
    spline interpolation)

**Discretization**
    Mesh based operators (e.g., interpolation, L2 projection)

**Benchmarks**
    Mesh and partitioning infrastructure of problems relevant to DTK

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
