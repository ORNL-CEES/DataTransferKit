ARG BASE=ubuntu:16.04
FROM $BASE

ARG NPROC=4

RUN apt-get update && apt-get install -y \
        build-essential \
        bc \
        gfortran \
        curl \
        git \
        wget \
        vim \
        emacs \
        autoconf \
        lcov \
        valgrind \
        ccache \
        cppcheck \
        libssl-dev \
        libpng-dev \
        libfreetype6-dev \
        libxft-dev \
        libsqlite3-dev \
        libbz2-dev \
        libcurl4-gnutls-dev \
        libatlas-base-dev \
        zlib1g-dev \
        python2.7-dev \
        ninja-build \
        doxygen \
        && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get install -y \
        gcc-7 \
        g++-7 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN curl -s https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py && rm get-pip.py && \
    pip install --upgrade --no-cache-dir pip && \
    pip install --no-cache-dir sphinx==1.4 sphinx_rtd_theme breathe==4.7

ENV PREFIX=/scratch
RUN mkdir -p ${PREFIX} && \
    cd ${PREFIX} && \
    mkdir archive && \
    mkdir source && \
    mkdir build && \
    mkdir install

# Install CMake
RUN export CMAKE_VERSION=3.13.3 && \
    export CMAKE_VERSION_SHORT=3.13 && \
    export CMAKE_URL=https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    export CMAKE_SCRIPT=cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    export CMAKE_PREFIX=/opt/cmake/${CMAKE_VERSION_SHORT} && \
    wget --quiet ${CMAKE_URL} --output-document=${CMAKE_SCRIPT} && \
    mkdir -p ${CMAKE_PREFIX} && \
    sh ${CMAKE_SCRIPT} --skip-license --prefix=${CMAKE_PREFIX} && \
    rm ${CMAKE_SCRIPT}
ENV PATH=/opt/cmake/3.13/bin:$PATH

# Install Clang/LLVM
RUN export LLVM_VERSION=6.0.1 && \
    export LLVM_URL=http://releases.llvm.org/${LLVM_VERSION}/clang+llvm-${LLVM_VERSION}-x86_64-linux-gnu-ubuntu-16.04.tar.xz && \
    export LLVM_ARCHIVE=clang+llvm-${LLVM_VERSION}-x86_64-linux-gnu-ubuntu-16.04.tar.xz && \
    export LLVM_PREFIX=/opt/llvm/6.0 && \
    wget --quiet ${LLVM_URL} --output-document=${LLVM_ARCHIVE} && \
    mkdir -p ${LLVM_PREFIX} && \
    tar -xvf ${LLVM_ARCHIVE} -C ${LLVM_PREFIX} --strip-components=1 && \
    echo "${LLVM_PREFIX}/lib" > /etc/ld.so.conf.d/llvm.conf && ldconfig && \
    rm -rf ${LLVM_ARCHIVE}
ENV PATH=/opt/llvm/6.0/bin:$PATH

# Install OpenMPI
RUN export OPENMPI_VERSION=2.1.6 && \
    export OPENMPI_VERSION_SHORT=2.1 && \
    export OPENMPI_SHA1=7a1d7f1b7efe2258bda3929b8b1729bfb7a51f08 && \
    export OPENMPI_URL=https://www.open-mpi.org/software/ompi/v${OPENMPI_VERSION_SHORT}/downloads/openmpi-${OPENMPI_VERSION}.tar.bz2 && \
    export OPENMPI_ARCHIVE=${PREFIX}/archive/openmpi-${OPENMPI_VERSION}.tar.bz2 && \
    export OPENMPI_SOURCE_DIR=${PREFIX}/source/openmpi/${OPENMPI_VERSION} && \
    export OPENMPI_BUILD_DIR=${PREFIX}/build/openmpi/${OPENMPI_VERSION} && \
    export OPENMPI_INSTALL_DIR=/opt/openmpi/${OPENMPI_VERSION_SHORT} && \
    wget --quiet ${OPENMPI_URL} --output-document=${OPENMPI_ARCHIVE} && \
    echo "${OPENMPI_SHA1} ${OPENMPI_ARCHIVE}" | sha1sum -c && \
    mkdir -p ${OPENMPI_SOURCE_DIR} && \
    tar -xf ${OPENMPI_ARCHIVE} -C ${OPENMPI_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${OPENMPI_BUILD_DIR} && \
    cd ${OPENMPI_BUILD_DIR} && \
    ${OPENMPI_SOURCE_DIR}/configure --prefix=${OPENMPI_INSTALL_DIR} && \
    make -j${NPROC} install && \
    rm -rf ${OPENMPI_ARCHIVE} && \
    rm -rf ${OPENMPI_BUILD_DIR} && \
    rm -rf ${OPENMPI_SOURCE_DIR}

ENV MPI_DIR=/opt/openmpi/2.1
# Put OPENMPI_DIR at the end of the path so that /ust/local/bin/mpiexec will
# overwrite it
ENV PATH=$PATH:${MPI_DIR}/bin


# install Boost
RUN export BOOST_VERSION=1.67.0 && \
    export BOOST_VERSION_UNDERSCORE=$(echo "$BOOST_VERSION" | sed -e "s/\./_/g") && \
    export BOOST_URL=https://dl.bintray.com/boostorg/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_UNDERSCORE}.tar.bz2 && \
    export BOOST_SHA256=2684c972994ee57fc5632e03bf044746f6eb45d4920c343937a465fd67a5adba && \
    export BOOST_ARCHIVE=${PREFIX}/archive/boost_${BOOST_VERSION_UNDERSCORE}.tar.bz2 && \
    export BOOST_SOURCE_DIR=${PREFIX}/source/boost/${BOOST_VERSION} && \
    export BOOST_BUILD_DIR=${PREFIX}/build/boost/${BOOST_VERSION} && \
    export BOOST_INSTALL_DIR=/opt/boost/${BOOST_VERSION} && \
    wget --quiet ${BOOST_URL} --output-document=${BOOST_ARCHIVE} && \
    echo "${BOOST_SHA256} ${BOOST_ARCHIVE}" | sha256sum -c && \
    mkdir -p ${BOOST_SOURCE_DIR} && \
    tar -xf ${BOOST_ARCHIVE} -C ${BOOST_SOURCE_DIR} --strip-components=1 && \
    cd ${BOOST_SOURCE_DIR} && \
    ./bootstrap.sh \
        --prefix=${BOOST_INSTALL_DIR} \
        && \
    echo "using mpi ;" >> project-config.jam && \
    ./b2 -j${NPROC} \
        --build-dir=${BOOST_BUILD_DIR} \
        hardcode-dll-paths=true dll-path=${BOOST_INSTALL_DIR}/lib \
        link=shared \
        variant=release \
        install \
        && \
    rm -rf ${BOOST_ARCHIVE} && \
    rm -rf ${BOOST_BUILD_DIR} && \
    rm -rf ${BOOST_SOURCE_DIR}

ENV BOOST_DIR=/opt/boost/1.67.0

# install HDF5
RUN export HDF5_VERSION=1.10.2 && \
    export HDF5_URL=http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.bz2 && \
    export HDF5_MD5=41fb9347801b546fba323523a1c1af51 && \
    export HDF5_ARCHIVE=${PREFIX}/archive/hdf5-${HDF5_VERSION}.tar.bz2 && \
    export HDF5_SOURCE_DIR=${PREFIX}/source/hdf5/${HDF5_VERSION} && \
    export HDF5_BUILD_DIR=${PREFIX}/build/hdf5/${HDF5_VERSION} && \
    export HDF5_INSTALL_DIR=/opt/hdf5/${HDF5_VERSION} && \
    wget --quiet ${HDF5_URL} --output-document=${HDF5_ARCHIVE} && \
    echo "${HDF5_MD5} ${HDF5_ARCHIVE}" | md5sum -c && \
    mkdir -p ${HDF5_SOURCE_DIR} && \
    tar -xf ${HDF5_ARCHIVE} -C ${HDF5_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${HDF5_BUILD_DIR} && \
    cd ${HDF5_BUILD_DIR} && \
    ${HDF5_SOURCE_DIR}/configure \
        --prefix=${HDF5_INSTALL_DIR} \
        --enable-shared \
        --disable-static \
        --enable-parallel \
        && \
    make -j${NPROC} install && \
    rm -rf ${HDF5_ARCHIVE} && \
    rm -rf ${HDF5_BUILD_DIR} && \
    rm -rf ${HDF5_SOURCE_DIR}

ENV HDF5_DIR=/opt/hdf5/1.10.2

# install NetCDF
RUN export NETCDF_VERSION=4.6.1 && \
    export NETCDF_URL=https://github.com/Unidata/netcdf-c/archive/v${NETCDF_VERSION}.tar.gz && \
    export NETCDF_ARCHIVE=${PREFIX}/archive/netcdf-${NETCDF_VERSION}.tar.gz && \
    export NETCDF_SOURCE_DIR=${PREFIX}/source/netcdf/${NETCDF_VERSION} && \
    export NETCDF_BUILD_DIR=${PREFIX}/build/netcdf/${NETCDF_VERSION} && \
    export NETCDF_INSTALL_DIR=/opt/netcdf/${NETCDF_VERSION} && \
    wget --quiet ${NETCDF_URL} --output-document=${NETCDF_ARCHIVE} && \
    mkdir -p ${NETCDF_SOURCE_DIR} && \
    tar -xf ${NETCDF_ARCHIVE} -C ${NETCDF_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${NETCDF_BUILD_DIR} && \
    cd ${NETCDF_BUILD_DIR} && \
    ${NETCDF_SOURCE_DIR}/configure \
        --prefix=${NETCDF_INSTALL_DIR} \
        --enable-netcdf-4 \
        --enable-shared \
        --disable-static \
        CC=${MPI_DIR}/bin/mpicc \
        CFLAGS="-I${HDF5_DIR}/include" \
        LDFLAGS="-L${HDF5_DIR}/lib -lhdf5" \
        && \
    make -j${NPROC} install && \
    rm -rf ${NETCDF_ARCHIVE} && \
    rm -rf ${NETCDF_BUILD_DIR} && \
    rm -rf ${NETCDF_SOURCE_DIR}

ENV NETCDF_DIR=/opt/netcdf/4.6.1


# download Trilinos
# Current hash corresponds to the merge of Kokkos 2.8.0
RUN export TRILINOS_HASH=2a24058eab4e932fe961717fa0dd860c2fcbb52b && \
    export TRILINOS_SHORT_HASH=2a24058e && \
    export TRILINOS_URL=https://github.com/trilinos/Trilinos/archive/${TRILINOS_HASH}.tar.gz && \
    export TRILINOS_ARCHIVE=${PREFIX}/archive/trilinos-${TRILINOS_HASH}.tar.gz && \
    export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_SHORT_HASH} && \
    export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_SHORT_HASH} && \
    export TRILINOS_INSTALL_DIR=${PREFIX}/install/trilinos/${TRILINOS_SHORT_HASH} && \
    wget --quiet ${TRILINOS_URL} --output-document=${TRILINOS_ARCHIVE} && \
    mkdir -p ${TRILINOS_SOURCE_DIR} && \
    tar -xf ${TRILINOS_ARCHIVE} -C ${TRILINOS_SOURCE_DIR} --strip-components=1 && \
    ln -s ${TRILINOS_SOURCE_DIR} ${PREFIX}/source/trilinos/release && \
    mkdir -p ${TRILINOS_BUILD_DIR} && \
    rm -rf ${TRILINOS_ARCHIVE}

ENV TRILINOS_DIR=/scratch/source/trilinos/release


# download Kokkos Profiling and Debugging Tools
# sh does not support arrays so we need to use bash
# Note: the commit hash provided below has been tagged as Version 2.7.24 on GitHub
RUN ["/bin/bash","-c","export KOKKOS_TOOLS_HASH=f5c30224dc620471713f1cd491931f93e8678cac && \
    export KOKKOS_TOOLS_SHORT_HASH=f5c3022 && \
    export KOKKOS_TOOLS_ARRAY=(\"kernel-filter\" \
                               \"kernel-logger\" \
                               \"memory-events\" \
                               \"memory-hwm-mpi\" \
                               \"memory-hwm\" \
                               \"memory-usage\" \
                               \"simple-kernel-timer-json\" \
                               \"simple-kernel-timer\" \
                               \"space-time-stack\" \
                              ) && \
    export KOKKOS_TOOLS_URL=https://github.com/kokkos/kokkos-tools/archive/${KOKKOS_TOOLS_HASH}.tar.gz && \
    export KOKKOS_TOOLS_ARCHIVE=${PREFIX}/archive/kokkos-tools-${KOKKOS_TOOLS_HASH}.tar.gz && \
    export KOKKOS_TOOLS_SOURCE_DIR=${PREFIX}/source/kokkos-tools/${KOKKOS_TOOLS_SHORT_HASH} && \
    wget --quiet ${KOKKOS_TOOLS_URL} --output-document=${KOKKOS_TOOLS_ARCHIVE} && \
    mkdir -p ${KOKKOS_TOOLS_SOURCE_DIR} && \
    tar -xf ${KOKKOS_TOOLS_ARCHIVE} -C ${KOKKOS_TOOLS_SOURCE_DIR} --strip-components=1 && \
    cd ${KOKKOS_TOOLS_SOURCE_DIR}/src/tools && \
    for tool in ${KOKKOS_TOOLS_ARRAY[@]}; do cd $tool && make && cd ../; done && \
    ln -s ${KOKKOS_TOOLS_SOURCE_DIR} ${PREFIX}/source/kokkos-tools/release && \
    rm -rf ${KOKKOS_TOOLS_ARCHIVE}"]

ENV KOKKOS_TOOLS_DIR=/scratch/source/kokkos-tools/release/src/tools

# Benchmark support library
ENV BENCHMARK_DIR=/opt/benchmark
RUN cd ${PREFIX} && \
    git clone https://github.com/google/benchmark.git -b v1.4.1 && \
    cd benchmark && \
    git clone https://github.com/google/googletest.git -b release-1.8.1 && \
    mkdir build && cd build && \
    cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=${BENCHMARK_DIR} .. && \
    make -j${NPROC} && make install && \
    cd ../.. && rm -rf benchmark

# append the option flag --allow-run-as-root to mpiexec
RUN echo '#!/usr/bin/env bash' > /usr/local/bin/mpiexec && \
    echo '${MPI_DIR}/bin/mpiexec --allow-run-as-root "$@"' >> /usr/local/bin/mpiexec && \
    chmod +x /usr/local/bin/mpiexec && \
    mpiexec --version

# setup vim
COPY .vimrc /root/.vimrc
# FIXME: workaround for CMake (< 3.11) to detect Boost.Python (>= 1.67) with version suffix
RUN ln -s ${BOOST_DIR}/lib/libboost_python27.so ${BOOST_DIR}/lib/libboost_python.so
RUN git clone https://github.com/VundleVim/Vundle.vim.git /root/.vim/bundle/Vundle.vim && \
    vim +PluginInstall +qall && \
    cd /root/.vim/bundle/YouCompleteMe && \
    BOOST_ROOT=${BOOST_DIR} ./install.py --clang-completer --system-libclang --system-boost
