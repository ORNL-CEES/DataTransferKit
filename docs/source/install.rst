Installation
============

This section provide guidelines for installing DataTransferKit and its TPLs.

Install third-party libraries
-----------------------------

The following third party libraries (TPLs) are used by DTK:

+------------------------+-------------------------------------+---------+
| Packages               | Dependency                          | Version |
+========================+=====================================+=========+
| Boost                  | Optional                            | 1.61.0  |
+------------------------+-------------------------------------+---------+
| BLAS/LAPACK            | Required                            | N/A     |
+------------------------+-------------------------------------+---------+
| MPI                    | Required                            | N/A     |
+------------------------+-------------------------------------+---------+
| Trilinos               | Required                            | 12.X    |
+------------------------+-------------------------------------+---------+
| Google Benchmark       | Required (if examples are enabled)  | 1.4     |
+------------------------+-------------------------------------+---------+

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
    $ spack install benchmark

Once installed, the module files for the packages must be loaded into the
environment by doing

.. code::

    $ spack load openblas
    $ spack load boost
    $ spack load openmpi
    $ spack load benchmark


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
        -D CMAKE_BUILD_TYPE=Release \
        -D TPL_ENABLE_MPI=ON \
        -D TPL_ENABLE_BLAS=ON \
        -D TPL_ENABLE_LAPACK=ON \
        -D TPL_ENABLE_Boost=ON \
        -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
        -D Tpetra_INST_INT_UNSIGNED_LONG=ON \
        -D Tpetra_INST_INT_LONG_LONG=OFF \
        -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
        -D Trilinos_EXTRA_REPOSITORIES="DataTransferKit" \
        -D Trilinos_ENABLE_DataTransferKit=ON \
        -D DataTransferKit_ENABLE_DBC=ON \
        -D DataTransferKit_ENABLE_TESTS=ON \
        -D DataTransferKit_ENABLE_EXAMPLES=ON \
        $EXTRA_ARGS \
        $TRILINOS_DIR

and run it from your build directory:

.. code::

    $ mkdir build && cd build
    $ ../do-configure

More install scripts can be found in ``scripts/``.

.. note::

    The above ``do-configure`` script may get outdated. You can always refer to
    ``scripts/docker_cmake`` which is used in the Jenkins CI builds and
    therefore is required to be always up-to-date.

Build this documentation
------------------------

Building documentation requires `sphinx <http://www.sphinx-doc.org>`_.
(Re)configure with ``-D DataTransferKit_ENABLE_ReadTheDocs=ON`` and run:

.. code::

    $ make docs

Open the ``index.html`` in the directory ``DataTransferKit/docs/html``.

Generate Doxygen documentation
------------------------------

Configure with ``-D DataTransferKit_ENABLE_Doxygen=ON`` and run:

.. code::

    $ make doxygen

Checkout ``DataTransferKit/docs/doxygen/html/index.html``.
