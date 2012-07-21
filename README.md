README
======

The DataTransferKit is a package providing parallel mesh search and
field transfer services.


DEPENDENCIES
------------

The DataTransferKit has direct dependecies on the following Trilinos packages:

* Teuchos
* Tpetra
* Shards
* Intrepid
* Zoltan

In addition, the mesh database MOAB is required. DTK has been tested
with a stable release version of MOAB 4.5. You can get MOAB
[here](http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB).

For parallel builds, an MPI implementation is also required. Both Open
MPI and MPICH have been tested.


CONFIGURE and BUILD
-------------------

The DataTransferKit uses the TriBITS cmake build system. Look
[here](https://github.com/CNERG/DataTransferKit/tree/dev/doc/build_notes)
for a sample configure script and additional configure/build
documentation.

The following compilers have been tested with DTK.

* gcc-4.7.0
* gcc-4.6.1


TEST
----

The ctest harness is used for testing. Tests for each subpackage can
be found in DataTransferKit/$SUBPACKAGE/test.


EXAMPLES
--------

Several examples are provided for using the DataTransferKit for
parallel search and transfer operations. See the example directory in
each subpackage.


DOCUMENTATION
-------------

User and developer documentation is provided by Doxygen
[here](http://cnerg.github.com/DataTransferKit/). A domain model
document can be found
[here](https://github.com/CNERG/DataTransferKit/tree/dev/doc/domain_model).

BUG REPORTING and ISSUE TRACKING
--------------------------------

Bug reporting and issue tracking are provided by GitHub. Please report
all bugs [here](https://github.com/CNERG/DataTransferKit/issues) by
creating a new issue.

AUTHORS
-------

The following people have made substantial contributions to the
development of DataTransferKit:

* [Stuart Slattery](http://github.com/sslattery)
