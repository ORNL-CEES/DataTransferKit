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
#include "DTK_DBC.hpp"

#include <Kokkos_Core.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <type_traits>

//---------------------------------------------------------------------------//
// TEST FUNCTORS
//---------------------------------------------------------------------------//
// Fill a view with an index.
template<class View>
class FillFunctor
{
  public:

    FillFunctor( View data )
        : _data( data )
    { /* ... */ }

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t i) const 
    { _data(i) = i; }

  private:

    View _data;
};

//---------------------------------------------------------------------------//
// Assign one 2D view to another.
template<class View1, class View2>
class AssignFunctor
{
  public:

    AssignFunctor( View1 view_1, View2 view_2 )
        : _view_1(view_1)
        , _view_2(view_2)
    {
        static_assert( 
            std::is_same<typename View1::traits::data_type,
                         typename View2::traits::data_type>::value,
            "View data types must be the same" );
        static_assert( 
            std::is_same<typename View1::traits::device_type,
                         typename View2::traits::device_type>::value,
            "View device types must be the same" );
        static_assert( static_cast<unsigned>(View1::Rank) == 2,
                       "View ranks must be 2" );
        static_assert( static_cast<unsigned>(View2::Rank) == 2,
                       "View ranks must be 2" );
        DTK_REQUIRE( view_1.extent(0) == view_2.extent(0) );
        DTK_REQUIRE( view_1.extent(1) == view_2.extent(1) );
        _extent = view_1.extent(1);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t i) const 
    { 
        for ( int n = 0; n < _extent; ++n ) 
        {
            _view_2(i,n) = _view_1(i,n);
        }
    }

  private:

    int _extent;
    View1 _view_1;
    View2 _view_2;
};

//---------------------------------------------------------------------------//
// Sum the values in a view.
template<class Scalar, class View>
class SumFunctor
{
  public:

    SumFunctor( View data )
        : _data( data )
    { /* ... */ }

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t i, Scalar& val) const 
    { val += _data(i); }

  private:

    View _data;
};

//---------------------------------------------------------------------------//
// TEST TEMPLATE DECLARATIONS
//---------------------------------------------------------------------------//
// Test creating a view and run a basic parallel for kernel.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, basic_for_kernel, Scalar, Device )
{
    // Create a view in the execution space.
    using ViewType = Kokkos::View<Scalar*,Device>;
    int size = 1000;
    ViewType data( "data", size );

    // Populate the view in the execution space.
    Kokkos::parallel_for( Kokkos::RangePolicy<typename Device::execution_space>(0,size),
                          FillFunctor<ViewType>(data) );

    // Mirror the view to the host space and check the results.
    typename ViewType::HostMirror host_data =
        Kokkos::create_mirror_view( data );
    for ( int i = 0; i < size; ++i )
    {
        TEST_EQUALITY( host_data(i), i );
    }
}

//---------------------------------------------------------------------------//
// Test assigning views with different layouts.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, layout_assign_kernel, Scalar, Device )
{
    // Create a view in the execution space.
    using ViewType1 = Kokkos::View<Scalar**,Kokkos::LayoutLeft,Device>;
    int size = 1000;
    int dim = 3;
    ViewType1 data_1( "data_1", size, dim );

    // Populate the first view.
    for ( int d = 0; d < dim; ++d )
    {
        auto sv = Kokkos::subview(data_1,Kokkos::ALL(),d);
        Kokkos::parallel_for( Kokkos::RangePolicy<typename Device::execution_space>(0,size), 
                              FillFunctor<decltype(sv)>(sv) );
    }

    // Create another view.
    using ViewType2 = Kokkos::View<Scalar**,Kokkos::LayoutRight,Device>;
    ViewType2 data_2( "data_2", size, dim );

    // Copy the first view into the second.
    Kokkos::parallel_for( 
        Kokkos::RangePolicy<typename Device::execution_space>(0,size),
        AssignFunctor<ViewType1,ViewType2>(data_1,data_2) );

    // Check the second view on the host.
    typename ViewType2::HostMirror host_data = 
        Kokkos::create_mirror_view( data_2 );
    for ( int i = 0; i < size; ++i )
    {
        for ( int d = 0; d < dim; ++d )
        {
            TEST_EQUALITY( host_data(i,d), i );
        }
    }
}

//---------------------------------------------------------------------------//
// Test creating a view and run a basic reduction for kernel.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( View, basic_reduce_kernel, Scalar, Device )
{
    // Create a view in the execution space.
    using ViewType = Kokkos::View<Scalar*,Device>;
    int size = 1000;
    ViewType data( "data", size );

    // Populate the view in the execution space.
    Kokkos::parallel_for( Kokkos::RangePolicy<typename Device::execution_space>(0,size),
                          FillFunctor<ViewType>(data) );

    // Sum the result.
    Scalar sum = Teuchos::ScalarTraits<Scalar>::zero();
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space>(0,size),
        SumFunctor<Scalar,ViewType>(data),
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
#define UNIT_TEST_GROUP( SCALAR, DEVICE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, basic_for_kernel, SCALAR, DEVICE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, layout_assign_kernel, SCALAR, DEVICE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( View, basic_reduce_kernel, SCALAR, DEVICE )    

//---------------------------------------------------------------------------//
// Create the unary function and instantiate the unit tests.
#define GROUP_INSTANT_1( DEVICE )                \
    UNIT_TEST_GROUP( double, DEVICE )
DTK_TEST_DEVICE_INSTANT( GROUP_INSTANT_1 )

#define GROUP_INSTANT_2( DEVICE )                \
    UNIT_TEST_GROUP( float, DEVICE )
DTK_TEST_DEVICE_INSTANT( GROUP_INSTANT_2 )

#define GROUP_INSTANT_3( DEVICE )                \
    UNIT_TEST_GROUP( int, DEVICE )
DTK_TEST_DEVICE_INSTANT( GROUP_INSTANT_3 )

//---------------------------------------------------------------------------//
// end tstKokkosView.cpp
//---------------------------------------------------------------------------//
