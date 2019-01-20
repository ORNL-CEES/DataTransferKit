/****************************************************************************
 * Copyright (c) 2012-2019 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
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
        fill2dKernel<Scalar>
            <<<num_blocks, block_size>>>( dtk_view.data(), dims[0], dims[1] );
    }
    else if ( 3 == num_dims )
    {
        fill3dKernel<Scalar>
            <<<num_blocks, block_size>>>( dtk_view, dims[0], dims[1], dims[2] );
    }

    // NOTE: Synchronization is required because kernels are launched
    // asynchronously. This is the equivalent of calling Kokkos::fence() after
    // a kernel. If these kernels were executed on a specific stream, only
    // that stream would need to be synchronized.
    cudaDeviceSynchronize();
}

//---------------------------------------------------------------------------//

} // namespace ViewTestCudaHelpers

#endif // end DTK_VIEWTESTCUDAHELPERS_DEF_HPP

//---------------------------------------------------------------------------//
// end ViewTestCudaHelpers_def.hpp
//---------------------------------------------------------------------------//
