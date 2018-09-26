from ubuntu:16.04

RUN apt-get update && apt-get install -y wget git doxygen python && \
    wget -q  https://bootstrap.pypa.io/get-pip.py -O get-pip.py && \
    python get-pip.py && rm get-pip.py && \
    pip install --upgrade --no-cache-dir pip && \
    pip install --no-cache-dir sphinx==1.7.5 sphinx_rtd_theme breathe==4.10 && \
    echo "deb http://apt.llvm.org/xenial/ llvm-toolchain-xenial-6.0 main" > /etc/apt/sources.list.d/llvm.list && \
    wget -q -O - http://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - && \
    apt-get update && apt-get install -y clang-format-6.0 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

