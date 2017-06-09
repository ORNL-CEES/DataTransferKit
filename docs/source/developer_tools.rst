Developer Tools
===============

Run DTK development environment in a Docker container
-----------------------------------------------------

To start a container from the DTK pre-built Docker image that is used in the
automated build on Jenkins, run:

.. code:: bash

    [host]$ cd docker
    [host]$ echo COMPOSE_PROJECT_NAME=$USER > .env # [optional] specify a project name
    [host]$ docker-compose pull # pull the most up-to-date version of the DTK base image
    [host]$ docker-compose up -d # start the container

This will mount the local DTK source directory into the container at
``${TRILINOS_DIR}/DataTransferKit``.  The environment variable ``TRILINOS_DIR``
is already defined and contains the path to a release version of Trilinos that
has been downloaded into the DTK base image.  We recommend you use a ``.env``
file to specify an alternate project name (the default being the directory name,
i.e. ``docker``).  This will let you run multiple isolated environments on a
single host.  Here the service name will be prefixed by your username which will
prevent interferences with other developers on the same system.

Then to launch an interactive Bash session inside that container, do:

.. code:: bash

    [host]$ docker-compose exec dtk_dev bash

Configure, build, and test as you would usually do:

.. code:: bash

    [container]$ cd $TRILINOS_DIR/DataTransferKit
    [container]$ mkdir build && cd build
    [container]$ ../scripts/docker_cmake
    [container]$ make -j<N>
    [container]$ ctest -j<N>

Do not forget to cleanup after yourself:

.. code:: bash

    [container]$ exit
    [host]$ docker-compose down # stop and remove the container


Make the container GPU-aware
----------------------------

Follow these instructions to launch containers leveraging the NVIDIA GPUs on the
host machine:

.. code:: bash

    [host]$ cd docker
    [host]$ nvidia/setup_nvidia_docker_compose.py # extend the regular Compose file to leverage GPUs
    [host]$ docker-compose build # add the CUDA development tools to the DTK base image
    [host]$ docker-compose up -d # as previously described

Do not forget to set the environment for CUDA before you configure:

.. code:: bash

    [container]$ # assuming you are in $TRILINOS_DIR/DataTransferKit/build
    [container]$ source ../scripts/set_kokkos_env.sh # set environment for CUDA
    [container]$ ../scripts/docker_cuda_cmake # configure
    [container]$ # now you may build and test


Code completion for Vim
-----------------------
Configure with ``-D DataTransferKit_ENABLE_YouCompleteMe=ON`` to generate a
``.ycm_extra_conf.py`` file at the root of your DTK source directory tree for
use with YouCompleteMe.
