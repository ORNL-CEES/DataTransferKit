#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=2}
# append the option flag --allow-run-as-root to mpiexec
mv /usr/bin/mpiexec /usr/bin/mpiexec.alias
echo '#!/usr/bin/env bash' > /usr/bin/mpiexc
echo 'mpiexec.alias --allow-run-as-root "$@"' >> /usr/bin/mpiexec
chmod +x /usr/bin/mpiexec
# cleanup workspace
rm -rf build
mkdir build
# reconfigure trilinos with dtk
TRILINOS_VERSION=12.4.2
export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_VERSION}
export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_VERSION}
export TRILINOS_INSTALL_DIR=${PREFIX}/install/trilinos/${TRILINOS_VERSION}

ln -sf ${WORKSPACE} ${TRILINOS_SOURCE_DIR}/DataTransferKit

cd ${TRILINOS_BUILD_DIR}
${TRILINOS_SOURCE_DIR}/DataTransferKit/docker/configure_trilinos.sh \
   -D Trilinos_EXTRA_REPOSITORIES:STRING=DataTransferKit \
   -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
   -D DataTransferKit_ENABLE_DBC:BOOL=ON \
   -D DataTransferKit_ENABLE_TESTS:BOOL=ON
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test

cp -r Testing ${WORKSPACE}/build/
