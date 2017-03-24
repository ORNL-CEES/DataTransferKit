#!/usr/bin/env bash

set -e

if [ -z "${CUDA_VERSION}" ]; then
    # Attempt to pull newer version of the image
    docker-compose pull ci
else
    # Update the configuration extension that makes CI container GPU-aware if
    # necessary.
    make setup_nvidia_docker
    # Attempt to pull a newer version of the image and rebuild (i.e. install
    # CUDA toolchain on top of it)
    docker-compose build --pull ci
fi

# View and validate the compose file
docker-compose config

# Request codecov to detect CI environment to pass through to docker
ci_env=`bash <(curl -s https://codecov.io/env)`
# Run the build within the container and remove it after it is done
# NOTE: passing BUILD_TAG from Jenkins as the project name.  This value will be
# prepended along with the 'ci' service name to the container on start up.
test ! -z $BUILD_TAG # exit with non zero status if env variable is not defined
docker-compose -p $BUILD_TAG run --rm $ci_env ci
