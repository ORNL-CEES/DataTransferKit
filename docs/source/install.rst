Installation
============

This section provide guidelines for installing DataTransferKit and its TPLs.

Install third-party libraries
-----------------------------

The following third party libraries (TPLs) are used by DTK:

+------------------------+------------+---------+
| Packages               | Dependency | Version |
+========================+============+=========+
| Boost                  | Optional   | 1.59.0  |
+------------------------+------------+---------+
| BLAS/LAPACK            | Required   | N/A     |
+------------------------+------------+---------+
| MPI                    | Optional   | N/A     |
+------------------------+------------+---------+
| Trilinos               | Required   | 12.X    |
+------------------------+------------+---------+

.. note::

    DTK is built as an external package of Trilinos. Thus, DTK needs the source of
    the Trilinos library rather than an installed version.

The dependencies of DataTransferKit may be built using `Spack
<https://github.com/llnl/spack>`_ package manager. You need to install the
following packages:

.. code::

    $ spack install openblas
    $ spack install boost
    $ spack install mpi

Once installed, the module files for the packages must be loaded into the
environment by doing

.. code::

    $ spack load openblas
    $ spack load boost
    $ spack load openmpi


DTKData repository
------------------

The DTKData repository contains mesh files used in DTK examples. To build the
examples and include the mesh files include the git submodule with your cloned
repository:

.. code::

    $ git submodule init
    $ git submodule update

Another way to achieve this is to pass the ``--recursive`` option to the ``git
clone`` command which will automatically initialize and update DTKData in the
DataTransferKit repository.

Building DTK
------------

DTK is configured and built using `TriBITS <https://tribits.org>`_.  DTK builds
within Trilinos effectively as an extension package.  First, link DTK into the
Trilinos main directory:

.. code::

    $ cd $TRILINOS_DIR
    $ ln -s $DTK_DIR

Create a ``do-configure`` script such as:

.. code-block:: bash

    EXTRA_ARGS=$@

    cmake \
        -D CMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR \
        -D CMAKE_BUILD_TYPE:STRING=DEBUG \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
        -D BUILD_SHARED_LIBS:BOOL=OFF \
        -D TPL_ENABLE_MPI:BOOL=ON \
        -D TPL_ENABLE_Boost:BOOL=ON \
        -D Boost_INCLUDE_DIRS:PATH=$BOOST_INCLUDE_DIR \
        -D TPL_ENABLE_Libmesh:BOOL=OFF \
        -D TPL_ENABLE_Netcdf:BOOL=OFF \
        -D TPL_ENABLE_MOAB:BOOL=OFF \
        -D TPL_ENABLE_BinUtils:BOOL=OFF \
        -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
        -D Trilinos_EXTRA_REPOSITORIES="DataTransferKit" \
        -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
        -D Trilinos_ENABLE_CXX11:BOOL=ON \
        -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
        -D DataTransferKit_ENABLE_DBC:BOOL=ON \
        -D DataTransferKit_ENABLE_TESTS:BOOL=ON \
        -D DataTransferKit_ENABLE_EXAMPLES:BOOL=ON \
        $EXTRA_ARGS \
        $TRILINOS_DIR

and run it from your build directory:

.. code::

    $ mkdir build && cd build
    $ ../do-configure

More install scripts can be found in ``scripts/``.

Build this documentation
------------------------

(Re)configure with ``-D DataTransferKit_ENABLE_ReadTheDocs=ON`` and run:

.. code::

    $ make docs

Open the ``index.html`` in the directory ``DataTransferKit/docs/html``.

Generate Doxygen documentation
------------------------------

Configure with ``-D DataTransferKit_ENABLE_Doxygen=ON`` and do:

.. code::

    $ make doxygen

Checkout ``DataTransferKit/docs/doxygen/html/index.html``.
