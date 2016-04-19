#!/usr/bin/env bash

TRILINOS_VERSION=12.4.2 && \
export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_VERSION} && \
export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_VERSION} && \
export TRILINOS_INSTALL_DIR=${PREFIX}/install/trilinos/${TRILINOS_VERSION} && \
useradd -m -s /bin/bash -N -u 1000 jovyan && \
chown jovyan ${PREFIX} -R && \
su jovyan <<EOF
ctest -j2 -VV -S ${PREFIX}/source/dtk/docker/TravisCI.cmake
EOF
