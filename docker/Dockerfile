ARG BASE=nvidia/cuda:10.2-devel
FROM $BASE

ARG NPROC=8

RUN apt-get update && apt-get install -y \
        build-essential \
        bc \
        gfortran \
        curl \
        git \
        wget \
        autoconf \
        lcov \
        ccache \
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
        doxygen \
        clang-format \
        clang \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PREFIX=/scratch
RUN mkdir -p ${PREFIX} && \
    cd ${PREFIX} && \
    mkdir archive && \
    mkdir source && \
    mkdir build && \
    mkdir install

# Install CMake
RUN export CMAKE_VERSION=3.16.4 && \
    export CMAKE_VERSION_SHORT=3.16 && \
    export CMAKE_URL=https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    export CMAKE_SCRIPT=cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    export CMAKE_PREFIX=/opt/cmake/${CMAKE_VERSION_SHORT} && \
    wget --quiet ${CMAKE_URL} --output-document=${CMAKE_SCRIPT} && \
    mkdir -p ${CMAKE_PREFIX} && \
    sh ${CMAKE_SCRIPT} --skip-license --prefix=${CMAKE_PREFIX} && \
    rm ${CMAKE_SCRIPT}
ENV PATH=/opt/cmake/3.16/bin:$PATH

# Install OpenMPI
RUN export OPENMPI_VERSION=4.0.2 && \
    export OPENMPI_VERSION_SHORT=4.0 && \
    export OPENMPI_SHA1=32ce3761288575fb8e4f6296c9105c3a25cf3235 && \
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

ENV MPI_DIR=/opt/openmpi/4.0
# Put OPENMPI_DIR at the end of the path so that /ust/local/bin/mpiexec will
# overwrite it
ENV PATH=$PATH:${MPI_DIR}/bin


# install Boost
RUN export BOOST_VERSION=1.67.0 && \
    export BOOST_VERSION_UNDERSCORE=$(echo "$BOOST_VERSION" | sed -e "s/\./_/g") && \
    export BOOST_URL=https://boostorg.jfrog.io/artifactory/main/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_UNDERSCORE}.tar.bz2 && \
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
RUN export HDF5_VERSION=1.10.6 && \
    export HDF5_URL=http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.bz2 && \
    export HDF5_MD5=03095102a6118c32a75a9b9b40be66f2 && \
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

ENV HDF5_DIR=/opt/hdf5/1.10.6

# install NetCDF
RUN export NETCDF_VERSION=4.7.3 && \
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

ENV NETCDF_DIR=/opt/netcdf/4.7.3

# download Trilinos 13.0.0
# Current hash has Kokkos 3.1
RUN export TRILINOS_HASH=9fec35276d846a667bc668ff4cbdfd8be0dfea08 && \
    export TRILINOS_SHORT_HASH=9fec3527 && \
    export TRILINOS_URL=https://github.com/trilinos/Trilinos/archive/${TRILINOS_HASH}.tar.gz && \
    export TRILINOS_ARCHIVE=${PREFIX}/archive/trilinos-${TRILINOS_HASH}.tar.gz && \
    export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_SHORT_HASH} && \
    export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_SHORT_HASH} && \
    export TRILINOS_INSTALL_DIR=/opt/trilinos && \
    export OMPI_CXX=${TRILINOS_SOURCE_DIR}/packages/kokkos/bin/nvcc_wrapper && \
    wget --quiet ${TRILINOS_URL} --output-document=${TRILINOS_ARCHIVE} && \
    mkdir -p ${TRILINOS_SOURCE_DIR} && \
    tar -xf ${TRILINOS_ARCHIVE} -C ${TRILINOS_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${TRILINOS_BUILD_DIR} && \
    cd ${TRILINOS_BUILD_DIR} && \
    cmake \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_STANDARD=14 \
      -D BUILD_SHARED_LIBS=ON \
      -D TPL_ENABLE_MPI=ON \
      -D TPL_ENABLE_BLAS=ON \
      -D TPL_ENABLE_LAPACK=ON \
      -D TPL_ENABLE_Boost=ON \
        -D Boost_INCLUDE_DIRS=${BOOST_DIR}/include \
        -D Boost_LIBRARY_DIRS=${BOOST_DIR}/lib \
      -D TPL_ENABLE_BoostLib=ON \
        -D BoostLib_INCLUDE_DIRS=${BOOST_DIR}/include \
        -D BoostLib_LIBRARY_DIRS=${BOOST_DIR}/lib \
      -D TPL_ENABLE_Netcdf=ON \
        -D Netcdf_INCLUDE_DIRS=$NETCDF_DIR/include \
        -D Netcdf_LIBRARY_DIRS=$NETCDF_DIR/lib \
      -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -D Trilinos_ENABLE_ALL_PACKAGES=OFF \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
      -D Trilinos_ENABLE_TESTS=OFF \
      -D Trilinos_ENABLE_EXAMPLES=OFF \
      -D Trilinos_ENABLE_OpenMP=ON \
      -D Trilinos_ENABLE_Teuchos=ON \
      -D Trilinos_ENABLE_Intrepid2=ON \
      -D Trilinos_ENABLE_Belos=ON \
      -D Trilinos_ENABLE_Stratimikos=ON \
      -D Trilinos_ENABLE_Thyra=ON \
      -D Trilinos_ENABLE_Tpetra=ON \
        -D Tpetra_INST_COMPLEX_DOUBLE=OFF \
        -D Tpetra_INST_COMPLEX_FLOAT=OFF \
        -D Tpetra_INST_SERIAL=ON \
        -D Tpetra_INST_OPENMP=ON \
      -D Kokkos_ENABLE_OpenMP=ON \
      -D TPL_ENABLE_CUDA=ON \
      -D Kokkos_ENABLE_Cuda=ON \
      -D Kokkos_ENABLE_Cuda_UVM=ON \
      -D Kokkos_ENABLE_Cuda_Lambda=ON \
      -D Kokkos_ARCH_VOLTA70=ON \
      -D Tpetra_INST_CUDA=ON \
      -D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR} \
    ${TRILINOS_SOURCE_DIR} && \
    make -j${NPROC} install && \
    export OMPI_CXX=${TRILINOS_INSTALL_DIR}/bin/nvcc_wrapper && \
    rm -rf ${TRILINOS_ARCHIVE} && \
    rm -rf ${TRILINOS_BUILD_DIR} && \
    rm -rf ${TRILINOS_SOURCE_DIR}

ENV TRILINOS_DIR=/opt/trilinos

# Benchmark support library
ENV BENCHMARK_DIR=/opt/benchmark
RUN cd ${PREFIX} && \
    git clone https://github.com/google/benchmark.git -b v1.5.0 && \
    cd benchmark && \
    git clone https://github.com/google/googletest.git -b release-1.10.0 && \
    mkdir build && cd build && \
    cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=${BENCHMARK_DIR} .. && \
    make -j${NPROC} && make install && \
    cd ../.. && rm -rf benchmark
