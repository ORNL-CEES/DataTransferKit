Installation
============

This section provide guidelines for installing DataTransferKit and its TPLs.

Install third-party libraries
-----------------------------

The following third party libraries (TPLs) are used by DTK:

+------------------------+------------+---------+
| Packages               | Dependency | Version |
+========================+============+=========+
| Trilinos               | Required   | 12.0    |
+------------------------+------------+---------+
| Boost                  | Optional   | 1.59.0  |
+------------------------+------------+---------+
| libMesh                | Optional   | 1.0     |
+------------------------+------------+---------+
| MOAB                   | Optional   | 4.9.1   |
+------------------------+------------+---------+
| BLAS/LAPACK            | Required   | N/A     |
+------------------------+------------+---------+

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

DTK is configured and built using TriBITS <https://tribits.org>. DTK builds
within Trilinos effectively as an extension package. First, link DTK into the Trilinos main directory:

.. code::

    $ cd $TRILINOS_DIR
    $ ln -s $DTK_DIR

DTK can then be configured along with the rest of Trilinos using
TriBITS. TriBITS is an extension of CMake which allows DTK to be configured
with scripts such as:

.. code-block:: bash

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
    -D Trilinos_ASSERT_MISSING_PACKAGES=OFF \
    -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D Trilinos_ENABLE_CXX11:BOOL=ON \
    -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
    -D Trilinos_ENABLE_DataTransferKitSTKMeshAdapters:BOOL=OFF \
    -D Trilinos_ENABLE_DataTransferKitMoabAdapters:BOOL=OFF \
    -D Trilinos_ENABLE_DataTransferKitLibmeshAdapters:BOOL=OFF \
    -D Trilinos_ENABLE_DataTransferKitClassicDTKAdapters:BOOL=ON \
    -D DataTransferKit_ENABLE_DBC:BOOL=ON \
    -D DataTransferKit_ENABLE_TESTS:BOOL=ON \
    -D DataTransferKit_ENABLE_EXAMPLES:BOOL=ON \
    $TRILINOS_DIR

More install scripts can be found in ``scripts/``.

Build this documentation
------------------------

(Re)configure with ``-D DataTransferKit_ENABLE_ReadTheDocs=ON`` and run:

.. code::

    $ make docs

Open the ``index.html`` in the directory ``docs/html``.



Developer Tools
===============

Run DTK development environment in a Docker container
-----------------------------------------------------

To start a container from the DTK pre-built Docker image that is used in the
automated build on Jenkins, run:

.. code:: bash

    [host]$ cd docker
    [host]$ docker-compose up -d

This will mount the local DTK source directory into the container at
``${TRILINOS_DIR}/DataTransferKit``.  The environment variable ``TRILINOS_DIR``
is already defined and contains the path to a release version of Trilinos that
has been downloaded into the DTK base image.

Then to launch an interactive Bash session inside that container, do:

.. code:: bash

    [host]$ docker exec -it dtk_dev bash

Configure, build, and test as you would usually do:

.. code::

    [container]$ cd $TRILINOS_DIR/DataTransferKit
    [container]$ mkdir build && cd build
    [container]$ ../scripts/docker_cmake
    [container]$ make -j<N>
    [container]$ ctest -j<N>

Do not forget to cleanup after yourself:

.. code:: bash

    [container]$ exit
    [host]$ docker-compose stop && docker-compose rm


Code completion for Vim
-----------------------
Configure with ``-D DataTransferKit_ENABLE_YouCompleteMe`` to generate a
``.ycm_extra_conf.py`` file at the root of your DTK source directory tree for
use with YouCompleteMe.
