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
 * \file   tstKokkosView.cpp
 * \author Stuart Slattery
 * \brief  Kokkos view unit tests.
 */
//---------------------------------------------------------------------------//

#include "DTK_UnitTestHelpers.hpp"

#include <Kokkos_Core.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>

//---------------------------------------------------------------------------//
// TEST TEMPLATE DECLARATIONS
//---------------------------------------------------------------------------//
// Test creating a view and run a basic parallel for kernel.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, basic_for_kernel, Scalar, ExecutionSpace )
{
    // Create a view in the execution space.
    using ViewType = Kokkos::View<Scalar*,ExecutionSpace>;
    int size = 1000;
    ViewType data( "data", size );

    // Populate the view in the execution space.
    Kokkos::parallel_for( size,
                          KOKKOS_LAMBDA(const size_t i){data(i) = i;} );

    // Mirror the view to the host space and check the results.
    typename ViewType::HostMirror host_data =
        Kokkos::create_mirror_view( data );
    for ( int i = 0; i < size; ++i )
    {
        TEST_EQUALITY( host_data(i), i );
    }
}

//---------------------------------------------------------------------------//
// Test creating a view and run a basic reduction for kernel.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, basic_reduce_kernel, Scalar, ExecutionSpace )
{
    // Create a view in the execution space.
    using ViewType = Kokkos::View<Scalar*,ExecutionSpace>;
    int size = 1000;
    ViewType data( "data", size );

    // Populate the view in the execution space.
    Kokkos::parallel_for( size,
                          KOKKOS_LAMBDA(const size_t i){data(i) = i;} );

    // Sum the result.
    Scalar sum = Teuchos::ScalarTraits<Scalar>::zero();
    Kokkos::parallel_reduce(
        size,
        KOKKOS_LAMBDA(const size_t i, Scalar& val){val += data(i);},
        sum );

    // Mirror the view to the host space and check the result.
    typename ViewType::HostMirror host_data =
        Kokkos::create_mirror_view( data );
    double test_sum = 0.0;
    for ( int i = 0; i < size; ++i )
    {
        test_sum += host_data(i);
    }
    TEST_EQUALITY( test_sum, sum );
}

//---------------------------------------------------------------------------//
// TEST TEMPLATE INSTANTIATIONS
//---------------------------------------------------------------------------//
// Create a unit test group.
#define UNIT_TEST_GROUP( SCALAR, SPACE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, basic_for_kernel, SCALAR, SPACE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, basic_reduce_kernel, SCALAR, SPACE )    

//---------------------------------------------------------------------------//
// Create the unary function and instantiate the unit tests.
#define GROUP_INSTANT_1( SPACE )                \
    UNIT_TEST_GROUP( double, SPACE )
DTK_TEST_SPACE_INSTANT( GROUP_INSTANT_1 )

#define GROUP_INSTANT_2( SPACE )                \
    UNIT_TEST_GROUP( float, SPACE )
DTK_TEST_SPACE_INSTANT( GROUP_INSTANT_2 )

#define GROUP_INSTANT_3( SPACE )                \
    UNIT_TEST_GROUP( int, SPACE )
DTK_TEST_SPACE_INSTANT( GROUP_INSTANT_3 )

//---------------------------------------------------------------------------//
// end tstKokkosView.cpp
//---------------------------------------------------------------------------//
