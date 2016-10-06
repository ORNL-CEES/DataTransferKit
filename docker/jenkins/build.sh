#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=8}
# append the option flag --allow-run-as-root to mpiexec
cat > /usr/local/bin/mpiexec <<\EOF
#!/usr/bin/env bash
/usr/bin/mpiexec --allow-run-as-root "$@"
EOF
chmod +x /usr/local/bin/mpiexec
# cleanup workspace
cd ${PREFIX}/source/dtk
rm -rf build
mkdir build
cd build
# configure trilinos with dtk
TRILINOS_VERSION=12.8.1
export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_VERSION}
#export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_VERSION}
export TRILINOS_INSTALL_DIR=${PREFIX}/install/trilinos/${TRILINOS_VERSION}
#cd ${TRILINOS_BUILD_DIR}
${TRILINOS_SOURCE_DIR}/DataTransferKit/docker/configure_trilinos.sh \
   -D Trilinos_EXTRA_REPOSITORIES:STRING=DataTransferKit \
   -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
   -D DataTransferKit_ENABLE_DBC:BOOL=ON \
   -D DataTransferKit_ENABLE_TESTS:BOOL=ON
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test

