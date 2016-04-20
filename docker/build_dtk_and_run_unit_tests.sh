#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=2}
# append the option flag --allow-run-as-root to mpiexec
mv ${MPI_DIR}/bin/mpiexec ${MPI_DIR}/bin/mpiexec.alias
echo '#!/usr/bin/env bash' > ${MPI_DIR}/bin/mpiexc
echo 'mpiexec.alias --allow-run-as-root "$@"' >> ${MPI_DIR}/bin/mpiexec
chmod +x ${MPI_DIR}/bin/mpiexec
# reconfigure trilinos with dtk
TRILINOS_VERSION=12.4.2
export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_VERSION}
export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_VERSION}
export TRILINOS_INSTALL_DIR=${PREFIX}/install/trilinos/${TRILINOS_VERSION}
cd ${TRILINOS_BUILD_DIR}
${PREFIX}/source/dtk/docker/configure_trilinos.sh \
   -D Trilinos_EXTRA_REPOSITORIES:STRING=DataTransferKit \
   -D Trilinos_ENABLE_DataTransferKit:BOOL=ON \
   -D DataTransferKit_ENABLE_DBC:BOOL=ON \
   -D DataTransferKit_ENABLE_TESTS:BOOL=ON
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} -V
