#!/user/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_DIR}/set_kokkos_env.sh

export CUDA_CXX_FLAGS="-g -arch=${GPU_ARCH} -lineinfo \
-Xcudafe --diag_suppress=conversion_function_not_usable \
-Xcudafe --diag_suppress=cc_clobber_ignored \
-Xcudafe --diag_suppress=code_is_unreachable"

export CUDA_EXTRA_ENABLE="\
-D TPL_ENABLE_CUDA=ON \
-D Kokkos_ENABLE_Cuda=ON \
-D Kokkos_ENABLE_Cuda_UVM=ON \
-D Tpetra_INST_CUDA=ON"
