/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_DETAILS_NEAREST_NEIGHBOR_OPERATOR_IMPL_HPP
#define DTK_DETAILS_NEAREST_NEIGHBOR_OPERATOR_IMPL_HPP

#include <DTK_DetailsDistributedSearchTreeImpl.hpp> // sendAcrossNetwork()
#include <DTK_DistributedSearchTree.hpp>

namespace DataTransferKit
{
namespace Details
{

template <typename DeviceType>
struct NearestNeighborOperatorImpl
{
    using ExecutionSpace = typename DeviceType::execution_space;

    static Kokkos::View<Nearest<DataTransferKit::Point> *, DeviceType>
    makeNearestNeighborQueries(
        Kokkos::View<Coordinate const **, DeviceType> target_points )
    {
        int const n_target_points = target_points.extent( 0 );
        Kokkos::View<Nearest<DataTransferKit::Point> *, DeviceType>
            nearest_queries( "nearest", n_target_points );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "setup_queries" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_target_points ),
            KOKKOS_LAMBDA( int i ) {
                nearest_queries( i ) = nearest(
                    Point{{target_points( i, 0 ), target_points( i, 1 ),
                           target_points( i, 2 )}} );
            } );
        Kokkos::fence();
        return nearest_queries;
    }

    template <typename View>
    static void
    pullSourceValues( MPI_Comm comm, View source_values,
                      Kokkos::View<int *, DeviceType> &buffer_indices,
                      Kokkos::View<int *, DeviceType> &buffer_ranks,
                      typename View::non_const_type &buffer_values )
    {
        static_assert(
            View::rank <= 2,
            "pullSourceValues() requires rank-1 or rank-2 view arguments" );
        int const n_exports = buffer_indices.extent( 0 );
        Distributor distributor( comm );
        int const n_imports = distributor.createFromSends( buffer_ranks );

        Kokkos::View<int *, DeviceType> export_target_indices( "target_indices",
                                                               n_exports );
        iota( export_target_indices );
        Kokkos::View<int *, DeviceType> import_target_indices( "target_indices",
                                                               n_imports );
        DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
            distributor, export_target_indices, import_target_indices );

        Kokkos::View<int *, DeviceType> export_source_indices = buffer_indices;
        Kokkos::View<int *, DeviceType> import_source_indices( "source_indices",
                                                               n_imports );
        DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
            distributor, export_source_indices, import_source_indices );

        Kokkos::View<int *, DeviceType> export_ranks( "ranks", n_exports );
        Kokkos::View<int *, DeviceType> import_ranks( "ranks", n_imports );
        int comm_rank;
        MPI_Comm_rank( comm, &comm_rank );
        Kokkos::deep_copy( export_ranks, comm_rank );
        DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
            distributor, export_ranks, import_ranks );

        buffer_indices = import_target_indices;
        buffer_ranks = import_ranks;
        Kokkos::realloc( buffer_values, n_imports,
                         source_values.dimension_1() );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "get_source_values" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
            KOKKOS_LAMBDA( int i ) {
                for ( int j = 0; j < (int)source_values.dimension_1(); ++j )
                    buffer_values( i, j ) =
                        source_values( import_source_indices( i ), j );
            } );
        Kokkos::fence();
    }

    template <typename View>
    static void
    pushTargetValues( MPI_Comm comm,
                      Kokkos::View<int *, DeviceType> const &buffer_indices,
                      Kokkos::View<int *, DeviceType> const &buffer_ranks,
                      View const &buffer_values, View target_values )
    {
        static_assert(
            View::rank <= 2,
            "pushTargetValues() requires rank-1 or rank-2 view arguments" );
        Distributor distributor( comm );
        int const n_imports = distributor.createFromSends( buffer_ranks );

        View export_source_values = buffer_values;
        View import_source_values( "source_values", n_imports,
                                   target_values.dimension_1() );
        DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
            distributor, export_source_values, import_source_values );

        Kokkos::View<int *, DeviceType> export_target_indices = buffer_indices;
        Kokkos::View<int *, DeviceType> import_target_indices( "target_indices",
                                                               n_imports );
        DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
            distributor, export_target_indices, import_target_indices );

        Kokkos::parallel_for(
            DTK_MARK_REGION( "set_target_values" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
            KOKKOS_LAMBDA( int i ) {
                for ( int j = 0; j < (int)target_values.dimension_1(); ++j )
                    target_values( import_target_indices( i ), j ) =
                        import_source_values( i, j );
            } );
        Kokkos::fence();
    }

    template <typename View>
    static typename View::non_const_type
    fetch( MPI_Comm comm, Kokkos::View<int const *, DeviceType> ranks,
           Kokkos::View<int const *, DeviceType> indices, View values )
    {
        DTK_REQUIRE( ranks.extent( 0 ) == indices.extent( 0 ) );

        Kokkos::View<int *, DeviceType> buffer_ranks =
            Kokkos::create_mirror( DeviceType(), ranks );
        Kokkos::deep_copy( buffer_ranks, ranks );

        Kokkos::View<int *, DeviceType> buffer_indices =
            Kokkos::create_mirror( DeviceType(), indices );
        Kokkos::deep_copy( buffer_indices, indices );

        typename View::non_const_type buffer_values( values.label() );

        pullSourceValues( comm, values, buffer_indices, buffer_ranks,
                          buffer_values );

        typename View::non_const_type values_out(
            values.label(), ranks.extent( 0 ), values.extent( 1 ) );

        pushTargetValues( comm, buffer_indices, buffer_ranks, buffer_values,
                          values_out );

        DTK_ENSURE( ( values_out.extent( 0 ) == ranks.extent( 0 ) ) &&
                    ( values_out.extent( 1 ) == values.extent( 1 ) ) );

        return values_out;
    }
};

} // namespace Details
} // namespace DataTransferKit

#endif
