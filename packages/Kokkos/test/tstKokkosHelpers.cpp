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
 * \file   tstKokkosHelpers.cpp
 * \author Stuart Slattery
 * \brief  Kokkos helpers unit tests.
 */
//---------------------------------------------------------------------------//

#include "DTK_DBC.hpp"
#include "DTK_KokkosHelpers.hpp"

#include <Kokkos_Core.hpp>

#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// TEST FUNCTORS
//---------------------------------------------------------------------------//
// Algorithm tags.
struct MinValue
{
};
struct MaxValue
{
};

//---------------------------------------------------------------------------//
// Fill function - specialize this for different fill operations.
template <class View, class Algorithm>
KOKKOS_INLINE_FUNCTION typename View::traits::value_type
fillFunction( Algorithm alg );

// Fill with the max of 0 and 1
template <class View>
KOKKOS_INLINE_FUNCTION
    typename View::traits::value_type fillFunction( MaxValue )
{
    typename View::traits::value_type zero = 0;
    typename View::traits::value_type one = 1;
    return DataTransferKit::KokkosHelpers::max( zero, one );
}

// Fill with the min of 0 and 1
template <class View>
KOKKOS_INLINE_FUNCTION
    typename View::traits::value_type fillFunction( MinValue )
{
    typename View::traits::value_type zero = 0;
    typename View::traits::value_type one = 1;
    return DataTransferKit::KokkosHelpers::min( zero, one );
}

//---------------------------------------------------------------------------//
// Fill a view with a given algorithm tag.
template <class View, class Algorithm>
class FillFunctor
{
  public:
    FillFunctor( View data )
        : _data( data )
    { /* ... */
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_t i ) const
    {
        _data( i ) = fillFunction<View>( Algorithm() );
    }

  private:
    View _data;
};

//---------------------------------------------------------------------------//
// TEST TEMPLATE DECLARATIONS
//---------------------------------------------------------------------------//
// Test helper functions.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( KokkosHelpers, helper_functions, Scalar,
                                   Node )
{
    // Get types.
    using DeviceType = typename Node::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;

    // Create a view in the execution space.
    using ViewType = Kokkos::View<Scalar *, DeviceType>;
    const int size = 10;
    ViewType data( "data", size );

    // Test the max
    {
        Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                              FillFunctor<ViewType, MaxValue>( data ) );
        Kokkos::fence();

        typename ViewType::HostMirror host_data =
            Kokkos::create_mirror_view( data );
        Kokkos::deep_copy( host_data, data );
        for ( int i = 0; i < size; ++i )
        {
            TEST_EQUALITY( host_data( i ), 1 );
        }
    }

    // Test the min
    {
        Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                              FillFunctor<ViewType, MinValue>( data ) );
        Kokkos::fence();

        typename ViewType::HostMirror host_data =
            Kokkos::create_mirror_view( data );
        Kokkos::deep_copy( host_data, data );
        for ( int i = 0; i < size; ++i )
        {
            TEST_EQUALITY( host_data( i ), 0 );
        }
    }
}

//---------------------------------------------------------------------------//
// TEST TEMPLATE INSTANTIATIONS
//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP_SN( SCALAR, NODE )                                     \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( KokkosHelpers, helper_functions,     \
                                          SCALAR, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SN( UNIT_TEST_GROUP_SN )

//---------------------------------------------------------------------------//
// end tstKokkosView.cpp
//---------------------------------------------------------------------------//
