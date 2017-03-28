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
 * \file   tstView.cpp
 * \author Stuart Slattery
 * \brief  DTK view unit tests.
 */
//---------------------------------------------------------------------------//

#include "DTK_View.hpp"

#if defined( KOKKOS_HAVE_CUDA )
#include "ViewTestCudaHelpers.hpp"
#endif // defined( KOKKOS_HAVE_CUDA )

#include <Kokkos_Core.hpp>

#include <Teuchos_UnitTestHarness.hpp>

#include <cassert>
#include <functional>
#include <type_traits>
#include <vector>

//---------------------------------------------------------------------------//
// View fill functions.
//---------------------------------------------------------------------------//
template <class Scalar, class ExectionSpace>
class FillView
{
  public:
    static void fill( DataTransferKit::View<Scalar> dtk_view,
                      const std::vector<unsigned> &dims );
};

//---------------------------------------------------------------------------//
// Serial implementation.
#if defined( KOKKOS_HAVE_SERIAL )
template <class Scalar>
class FillView<Scalar, Kokkos::Serial>
{
  public:
    static void fill( DataTransferKit::View<Scalar> dtk_view,
                      const std::vector<unsigned> &dims )
    {
        int num_dims = dims.size();

        if ( 1 == num_dims )
        {
            for ( unsigned i = 0; i < dims[0]; ++i )
            {
                dtk_view[i] = i;
            }
        }
        else if ( 2 == num_dims )
        {
            for ( unsigned i = 0; i < dims[0]; ++i )
            {
                for ( unsigned j = 0; j < dims[1]; ++j )
                {
                    dtk_view[j * dims[0] + i] = i + j;
                }
            }
        }
        else if ( 3 == num_dims )
        {
            for ( unsigned i = 0; i < dims[0]; ++i )
            {
                for ( unsigned j = 0; j < dims[1]; ++j )
                {
                    for ( unsigned k = 0; k < dims[2]; ++k )
                    {
                        dtk_view[k * dims[0] * dims[1] + j * dims[0] + i] =
                            i + j + k;
                    }
                }
            }
        }
    }
};
#endif // defined( KOKKOS_HAVE_SERIAL )

//---------------------------------------------------------------------------//
// OpenMP implementation.
#if defined( KOKKOS_HAVE_OPENMP )
template <class Scalar>
class FillView<Scalar, Kokkos::OpenMP>
{
  public:
    static void fill( DataTransferKit::View<Scalar> dtk_view,
                      const std::vector<unsigned> &dims )
    {
        int num_dims = dims.size();

        if ( 1 == num_dims )
        {
#pragma omp parallel for
            for ( unsigned i = 0; i < dims[0]; ++i )
            {
                dtk_view[i] = i;
            }
        }
        else if ( 2 == num_dims )
        {
#pragma omp parallel for collapse( 2 )
            for ( unsigned i = 0; i < dims[0]; ++i )
            {
                for ( unsigned j = 0; j < dims[1]; ++j )
                {
                    dtk_view[j * dims[0] + i] = i + j;
                }
            }
        }
        else if ( 3 == num_dims )
        {
#pragma omp parallel for collapse( 3 )
            for ( unsigned i = 0; i < dims[0]; ++i )
            {
                for ( unsigned j = 0; j < dims[1]; ++j )
                {
                    for ( unsigned k = 0; k < dims[2]; ++k )
                    {
                        dtk_view[k * dims[0] * dims[1] + j * dims[0] + i] =
                            i + j + k;
                    }
                }
            }
        }
    }
};
#endif // defined ( KOKKOS_HAVE_OPENMP )

//---------------------------------------------------------------------------//
// Cuda implementation.
#if defined( KOKKOS_HAVE_CUDA )
template <class Scalar>
class FillView<Scalar, Kokkos::Cuda>
{
  public:
    static void fill( DataTransferKit::View<Scalar> dtk_view,
                      const std::vector<unsigned> &dims )
    {
        ViewTestCudaHelpers::fillViewCuda( dtk_view, dims );
    }
};
#endif // defined( KOKKOS_HAVE_CUDA )

//---------------------------------------------------------------------------//
// TEST TEMPLATE DECLARATIONS
//---------------------------------------------------------------------------//
// Test creating a 1d view and run a basic parallel for kernel.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, 1d_view, Scalar, Node )
{
    // Get types.
    using DeviceType = typename Node::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;

    // Create a 1d view in the execution space.
    using ViewType = Kokkos::View<Scalar *, Kokkos::LayoutLeft, ExecutionSpace>;
    std::vector<unsigned> dims = {9};
    ViewType data( "1d_data", dims[0] );

    // Create a flat DTK view.
    DataTransferKit::View<Scalar> dtk_view( data );

    // Fill it with data.
    FillView<Scalar, ExecutionSpace>::fill( dtk_view, dims );

    // Mirror the view to the host space.
    typename ViewType::HostMirror host_data =
        Kokkos::create_mirror_view( data );
    Kokkos::deep_copy( host_data, data );

    // Check the results.
    for ( unsigned i = 0; i < dims[0]; ++i )
    {
        TEST_EQUALITY( host_data( i ), i );
    }
}

//---------------------------------------------------------------------------//
// Test creating a 1d view and run a basic parallel for kernel.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, 2d_view, Scalar, Node )
{
    // Get types.
    using DeviceType = typename Node::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;

    // Create a 2d view in the execution space.
    using ViewType =
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, ExecutionSpace>;
    std::vector<unsigned> dims = {9, 5};
    ViewType data( "2d_data", dims[0], dims[1] );

    // Create a flat DTK view.
    DataTransferKit::View<Scalar> dtk_view( data );

    // Fill it with data.
    FillView<Scalar, ExecutionSpace>::fill( dtk_view, dims );

    // Mirror the view to the host space.
    typename ViewType::HostMirror host_data =
        Kokkos::create_mirror_view( data );
    Kokkos::deep_copy( host_data, data );

    // Check the results.
    for ( unsigned j = 0; j < dims[1]; ++j )
    {
        for ( unsigned i = 0; i < dims[0]; ++i )
        {
            TEST_EQUALITY( host_data( i, j ), i + j );
        }
    }
}

//---------------------------------------------------------------------------//
// Test creating a 1d view and run a basic parallel for kernel.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, 3d_view, Scalar, Node )
{
    // Get types.
    using DeviceType = typename Node::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;

    // Create a 3d view in the execution space.
    using ViewType =
        Kokkos::View<Scalar ***, Kokkos::LayoutLeft, ExecutionSpace>;
    std::vector<unsigned> dims = {9, 5, 2};
    ViewType data( "3d_data", dims[0], dims[1], dims[2] );

    // Create a flat DTK view.
    DataTransferKit::View<Scalar> dtk_view( data );

    // Fill it with data.
    FillView<Scalar, ExecutionSpace>::fill( dtk_view, dims );

    // Mirror the view to the host space.
    typename ViewType::HostMirror host_data =
        Kokkos::create_mirror_view( data );
    Kokkos::deep_copy( host_data, data );

    // Check the results.
    for ( unsigned k = 0; k < dims[2]; ++k )
    {
        for ( unsigned j = 0; j < dims[1]; ++j )
        {
            for ( unsigned i = 0; i < dims[0]; ++i )
            {
                TEST_EQUALITY( host_data( i, j, k ), i + j + k );
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Test creating an empty view and call the default constructor.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, empty_view, Scalar, Node )
{
    DataTransferKit::View<Scalar> dtk_view;
    TEST_EQUALITY( dtk_view.size(), 0 );
    TEST_ASSERT( dtk_view.data() == nullptr );
}

//---------------------------------------------------------------------------//
// TEST TEMPLATE INSTANTIATIONS
//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP_SN( SCALAR, NODE )                                     \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, 1d_view, SCALAR, NODE )        \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, 2d_view, SCALAR, NODE )        \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, 3d_view, SCALAR, NODE )        \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, empty_view, SCALAR, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SN( UNIT_TEST_GROUP_SN )

//---------------------------------------------------------------------------//
// end tstView.cpp
//---------------------------------------------------------------------------//
