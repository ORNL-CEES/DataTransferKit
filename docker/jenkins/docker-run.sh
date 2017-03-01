#!/usr/bin/env bash

set -e

# Attempt to pull newer version of the image
docker-compose pull
# Build image if ci service specified (e.g. when extending the base image for
# NVIDIA Docker plugin)
docker-compose build --pull

# Request codecov to detect CI environment to pass through to docker
ci_env=`bash <(curl -s https://codecov.io/env)`
# Run the build within the container and remove it after it is done
# NOTE: passing BUILD_TAG from Jenkins as the project name.  This value will be
# prepended along with the 'ci' service name to the container on start up.
test ! -z $BUILD_TAG # exit with non zero status if env variable is not defined
docker-compose -p $BUILD_TAG run --rm $ci_env ci
