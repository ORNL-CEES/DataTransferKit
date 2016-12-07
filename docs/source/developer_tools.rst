Developer Tools
===============

Run DTK development environment in a Docker container
-----------------------------------------------------

To start a container from the DTK pre-built Docker image that is used in the
automated build on Jenkins, run:

.. code:: bash

    [host]$ cd docker
    [host]$ docker-compose pull # pull the most up-to-date version of the DTK base image
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
