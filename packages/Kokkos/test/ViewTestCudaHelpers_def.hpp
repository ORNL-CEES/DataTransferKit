//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file   ViewTestCudaHelpers_def.hpp
 * \author Stuart Slattery
 * \brief  tstView.cpp cuda test helpers
 */
//---------------------------------------------------------------------------//

#ifndef DTK_VIEWTESTCUDAHELPERS_DEF_HPP
#define DTK_VIEWTESTCUDAHELPERS_DEF_HPP

#include <cuda_runtime.h>

namespace ViewTestCudaHelpers
{

//---------------------------------------------------------------------------//
// 1d fill kernel.
template <class Scalar>
__global__ void fill1dKernel( DataTransferKit::View<Scalar> view,
                              const int dim1 )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < dim1 )
    {
        assert( idx < view.size() );
        view[idx] = idx;
    }
}

//---------------------------------------------------------------------------//
// 2d fill kernel
template <class Scalar>
__global__ void fill2dKernel( Scalar *data, const int dim1, const int dim2 )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < dim1 )
    {
        for ( int j = 0; j < dim2; ++j )
        {
            int index = j * dim1 + idx;
            assert( index < dim1 * dim2 );
            data[index] = idx + j;
        }
    }
}

//---------------------------------------------------------------------------//
// 3d fill kernel
template <class Scalar>
__global__ void fill3dKernel( DataTransferKit::View<Scalar> view,
                              const int dim1, const int dim2, const int dim3 )
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if ( idx < dim1 )
    {
        for ( int j = 0; j < dim2; ++j )
        {
            for ( int k = 0; k < dim3; ++k )
            {
                int index = k * dim1 * dim2 + j * dim1 + idx;
                assert( index < view.size() );
                view[index] = idx + j + k;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Fill function.
template <class Scalar>
void fillViewCuda( DataTransferKit::View<Scalar> dtk_view,
                   const std::vector<unsigned> &dims )
{
    int num_dims = dims.size();

    int block_size = 64;
    int num_blocks = dims[0] / block_size;
    if ( dims[0] % block_size > 0 )
        ++num_blocks;

    if ( 1 == num_dims )
    {
        fill1dKernel<Scalar><<<num_blocks, block_size>>>( dtk_view, dims[0] );
    }
    else if ( 2 == num_dims )
    {
        fill2dKernel<Scalar><<<num_blocks, block_size>>>( dtk_view.data(),
                                                          dims[0], dims[1] );
    }
    else if ( 3 == num_dims )
    {
        fill3dKernel<Scalar><<<num_blocks, block_size>>>( dtk_view, dims[0],
                                                          dims[1], dims[2] );
    }

    // NOTE: Synchronization is required because kernels are launched
    // asynchronously. This is the equivalent of calling Kokkos::fence() after
    // a kernel. If these kernels were executed on a specific stream, only
    // that stream would need to be synchronized.
    cudaDeviceSynchronize();
}

//---------------------------------------------------------------------------//

} // end namespace ViewTestCudaHelpers

#endif // end DTK_VIEWTESTCUDAHELPERS_DEF_HPP

//---------------------------------------------------------------------------//
// end ViewTestCudaHelpers_def.hpp
//---------------------------------------------------------------------------//
