#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=2}
# append the option flag --allow-run-as-root to mpiexec
cat > /usr/local/bin/mpiexec <<\EOF
#!/usr/bin/env bash
/usr/bin/mpiexec --allow-run-as-root "$@"
EOF
chmod +x /usr/local/bin/mpiexec
# configure trilinos with dtk
TRILINOS_VERSION=12.8.1
export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_VERSION}
export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_VERSION}
export TRILINOS_INSTALL_DIR=${PREFIX}/install/trilinos/${TRILINOS_VERSION}
cd ${TRILINOS_BUILD_DIR}
${PREFIX}/source/dtk/docker/configure_trilinos.sh
# build
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} -V
