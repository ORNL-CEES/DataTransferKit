/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <ArborX.hpp>
#include <DTK_DetailsNearestNeighborOperatorImpl.hpp> // fetch

#include <Teuchos_Array.hpp>
#include <Teuchos_UnitTestHarness.hpp>

template <
    typename View,
    typename std::enable_if<
        std::is_same<
            typename View::traits::memory_space,
            typename View::traits::host_mirror_space::memory_space>::value,
        int>::type = 0>
Teuchos::ArrayView<typename View::const_value_type> toArray( View const &v )
{
    return Teuchos::ArrayView<typename View::const_value_type>( v.data(),
                                                                v.size() );
}

template <
    typename View,
    typename std::enable_if<
        !std::is_same<
            typename View::traits::memory_space,
            typename View::traits::host_mirror_space::memory_space>::value,
        int>::type = 0>
Teuchos::Array<typename View::non_const_value_type> toArray( View const &v )
{
    auto v_h = Kokkos::create_mirror_view( v );
    Kokkos::deep_copy( v_h, v );

    return Teuchos::Array<typename View::non_const_value_type>(
        Teuchos::ArrayView<typename View::const_value_type>( v_h.data(),
                                                             v_h.size() ) );
}

// NOTE Shameless hack to get DeviceType from test driver because I don't want
// to deduce it
template <typename DeviceType>
struct Helper
{
    template <typename View1, typename View2>
    static void checkSendAcrossNetwork( MPI_Comm comm, View1 const &ranks,
                                        View2 const &v_exp, View2 const &v_ref,
                                        bool &success,
                                        Teuchos::FancyOStream &out )
    {
        ArborX::Details::Distributor<DeviceType> distributor( comm );
        distributor.createFromSends( typename DeviceType::execution_space{},
                                     ranks );

        // NOTE here we assume that the reference solution is sized properly
        auto v_imp =
            Kokkos::create_mirror( typename View2::memory_space(), v_ref );

        ArborX::Details::DistributedTreeImpl<DeviceType>::sendAcrossNetwork(
            typename DeviceType::execution_space{}, distributor, v_exp, v_imp );

        // FIXME not sure why I need that guy but I do get a bus error when it
        // is not here...
        Kokkos::fence();

        TEST_COMPARE_ARRAYS( toArray( v_imp ), toArray( v_ref ) );
    }

    template <typename View1, typename View2>
    static void checkFetch( MPI_Comm comm, View1 const &ranks,
                            View1 const &indices, View2 const &v_exp,
                            View2 const &v_ref, bool &success,
                            Teuchos::FancyOStream &out )
    {
        auto v_imp = DataTransferKit::Details::NearestNeighborOperatorImpl<
            DeviceType>::fetch( comm, ranks, indices, v_exp );

        TEST_COMPARE_ARRAYS( toArray( v_imp ), toArray( v_ref ) );
    }
};

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DetailsDistributedTreeImpl,
                                   send_across_network, DeviceType )
{
    using ExecutionSpace = typename DeviceType::execution_space;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    int comm_size;
    MPI_Comm_size( comm, &comm_size );

    int const DIM = 3;

    // send 1 packet to rank k
    // receive comm_size packets
    Kokkos::View<int **, DeviceType> u_exp( "u_exp", comm_size, DIM );
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size ),
                          KOKKOS_LAMBDA( int i ) {
                              for ( int j = 0; j < DIM; ++j )
                                  u_exp( i, j ) = i + j * comm_rank;
                          } );
    Kokkos::fence();

    Kokkos::View<int *, DeviceType> ranks_u( "", comm_size );
    ArborX::iota( ExecutionSpace{}, ranks_u, 0 );

    Kokkos::View<int **, DeviceType> u_ref( "u_ref", comm_size, DIM );
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size ),
                          KOKKOS_LAMBDA( int i ) {
                              for ( int j = 0; j < DIM; ++j )
                                  u_ref( i, j ) = comm_rank + i * j;
                          } );
    Kokkos::fence();

    Helper<DeviceType>::checkSendAcrossNetwork( comm, ranks_u, u_exp, u_ref,
                                                success, out );

    // send k packets to rank k
    // receive k*comm_size packets
    Kokkos::View<int *, DeviceType> tn( "tn", comm_size + 1 );
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size + 1 ),
        KOKKOS_LAMBDA( int i ) { tn( i ) = ( i * ( i - 1 ) ) / 2; } );
    Kokkos::fence();

    Kokkos::View<int **, DeviceType> v_exp( "v_exp", ArborX::lastElement( tn ),
                                            DIM );
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size ),
                          KOKKOS_LAMBDA( int i ) {
                              for ( int j = tn( i ); j < tn( i + 1 ); ++j )
                                  for ( int k = 0; k < DIM; ++k )
                                      v_exp( j, k ) = i * k;
                          } );
    Kokkos::fence();

    Kokkos::View<int *, DeviceType> ranks_v( "", ArborX::lastElement( tn ) );
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size ),
                          KOKKOS_LAMBDA( int i ) {
                              for ( int j = tn( i ); j < tn( i + 1 ); ++j )
                                  ranks_v( j ) = i;
                          } );
    Kokkos::fence();

    Kokkos::View<int **, DeviceType> v_ref( "v_ref", comm_size * comm_rank,
                                            DIM );
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size * comm_rank ),
        KOKKOS_LAMBDA( int i ) {
            for ( int j = 0; j < DIM; ++j )
                v_ref( i, j ) = j * comm_rank;
        } );
    Kokkos::fence();

    Helper<DeviceType>::checkSendAcrossNetwork( comm, ranks_v, v_exp, v_ref,
                                                success, out );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DetailsNearestNeighborOperatorImpl, fetch,
                                   DeviceType )
{
    using ExecutionSpace = typename DeviceType::execution_space;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    int comm_size;
    MPI_Comm_size( comm, &comm_size );

    // make communicaton plan
    Kokkos::View<int *, DeviceType> indices( "indices", comm_size );
    Kokkos::View<int *, DeviceType> ranks( "ranks", comm_size );
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size ),
                          KOKKOS_LAMBDA( int i ) {
                              indices( i ) = ( comm_rank + i ) % comm_size;
                              ranks( i ) = comm_size - 1 - comm_rank;
                          } );
    Kokkos::fence();

    // v(i) <-- k*comm_size+i (index i, rank k)
    Kokkos::View<int *, DeviceType> v_exp( "v", comm_size );
    ArborX::iota( ExecutionSpace{}, v_exp, comm_rank * comm_size );

    Kokkos::View<int *, DeviceType> v_ref( "v_ref", comm_size );
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size ),
                          KOKKOS_LAMBDA( int i ) {
                              v_ref( i ) =
                                  ranks( i ) * comm_size + indices( i );
                          } );
    Kokkos::fence();

    Helper<DeviceType>::checkFetch( comm, ranks, indices, v_exp, v_ref, success,
                                    out );

    // w(i, j) <-- k*comm_size*DIM+i+j*comm_size (index i, index j, rank k)
    int const DIM = 2;
    Kokkos::View<int **, DeviceType> w_exp( "w", comm_size, DIM );
    for ( int i = 0; i < DIM; ++i )
        ArborX::iota( ExecutionSpace{},
                      Kokkos::subview( w_exp, Kokkos::ALL, i ),
                      i * comm_size + comm_rank * comm_size * DIM );

    Kokkos::View<int **, DeviceType> w_ref( "w_ref", comm_size, DIM );
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, comm_size ),
                          KOKKOS_LAMBDA( int i ) {
                              for ( int j = 0; j < DIM; ++j )
                                  w_ref( i, j ) = ranks( i ) * comm_size * DIM +
                                                  indices( i ) + j * comm_size;
                          } );
    Kokkos::fence();

    Helper<DeviceType>::checkFetch( comm, ranks, indices, w_exp, w_ref, success,
                                    out );
}

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        DetailsDistributedTreeImpl, send_across_network, DeviceType##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DetailsNearestNeighborOperatorImpl,  \
                                          fetch, DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
