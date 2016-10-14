Installation
============

This section provide guidelines for installing DataTransferKit and its TPLs.

Install third-party libraries
-----------------------------

+------------------------+------------+---------+
| Packages               | Dependency | Version |
+========================+============+=========+
| Boost                  | Required   | 1.59.0? |
+------------------------+------------+---------+
| Trilinos               | Required   | 11.14   |
+------------------------+------------+---------+

Just a place holder. You paste and copy whatever you want here.


DTKData repository
------------------
Instead of creating a symlink to DTKData at the root of the DataTransferKit
source directory, you may use:

.. code::

    $ git submodule init
    $ git submodule update

Another way to achieve this is to pass the ``--recursive`` option to the ``git
clone`` command which will automatically initialize and update DTKData in the
DataTransferKit repository.


Build this documentation
------------------------

(Re)configure with ``-D DataTransferKit_ENABLE_ReadTheDocs=ON``

and try:

.. code::

    make docs

Open the ``index.html`` in the directory ``docs/html``.



Developer's corner
==================

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
