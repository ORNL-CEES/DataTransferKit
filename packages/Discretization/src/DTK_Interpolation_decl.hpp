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

#ifndef DTK_INTERPOLATION_DECL_HPP
#define DTK_INTERPOLATION_DECL_HPP

#include "DTK_ConfigDefs.hpp"
#include <ArborX.hpp>
#include <DTK_FE.hpp>
#include <DTK_FETypes.h>
#include <DTK_InterpolationFunctor.hpp>
#include <DTK_Mesh.hpp>
#include <DTK_PointSearch.hpp>
#include <DTK_Topology.hpp>

#include <Intrepid2_FunctionSpaceTools.hpp>

#include <mpi.h>

#include <array>
#include <string>

namespace DataTransferKit
{
/**
 * This class performs an interpolation for a set of given points in a given
 * mesh.
 */
template <typename DeviceType>
class Interpolation
{
  public:
    /**
     * Constructor.
     * @param comm
     * @param mesh mesh of the domain of interest
     * @param points_coordinates coordinates in the physical frame of the points
     * that we are looking for (n phys points, dim)
     * @param cell_dof_ids degrees of freedom indices associated to each cell (n
     * cells * n dofs per cell)
     * @param fe_type type of the finite element (DTK_HGRAD, DTK_HDIV, or
     * DTK_CURL)
     */
    Interpolation( MPI_Comm comm, Mesh<DeviceType> const &mesh,
                   Kokkos::View<Coordinate **, DeviceType> points_coordinates,
                   Kokkos::View<LocalOrdinal *, DeviceType> cell_dof_ids,
                   DTK_FEType fe_type );

    /**
     * This function performs the interpolation.
     * @param [in] X (n dofs, n fields)
     * @param [out] Y (n phys points, n fields)
     * @return View of size Y.extent(0) with the ID associated associated to
     * each physical points. This can be used to know if a point was not found
     * and which one it was.
     */
    template <typename Scalar>
    Kokkos::View<int *, DeviceType>
    apply( Kokkos::View<Scalar **, DeviceType> X,
           Kokkos::View<Scalar **, DeviceType> Y );

  private:
    void filter_dofs_ids(
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
        Kokkos::View<LocalOrdinal *, DeviceType> cell_dof_ids,
        DTK_FEType fe_type );

    /**
     * Helper function that calls Functor::Interpolation.
     */
    template <typename Scalar, typename FEOpType>
    void interpolate( Kokkos::View<Coordinate **, DeviceType> ref_points,
                      Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids,
                      Kokkos::View<Scalar **, DeviceType> X,
                      Kokkos::View<Scalar **, DeviceType> Y );

    /**
     * Helper function that calls Functor::HgradInterpolation.
     */
    template <typename Scalar, typename FEOpType>
    void
    hgradInterpolate( Kokkos::View<Coordinate **, DeviceType> ref_points,
                      Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids,
                      Kokkos::View<Scalar **, DeviceType> X,
                      Kokkos::View<Scalar **, DeviceType> Y );

    template <typename Scalar>
    void interpolateDispatch( FE fe, unsigned int fe_id,
                              Kokkos::View<Scalar **, DeviceType> X,
                              Kokkos::View<Scalar **, DeviceType> Y_fe );

    PointSearch<DeviceType> _point_search;

    /**
     * Dofs ids associated to each node of the cells.
     */
    std::array<Kokkos::View<LocalOrdinal **, DeviceType>, DTK_N_TOPO> _dofs_ids;

    /**
     * Map between the finite element index and the finite element basis.
     */
    std::array<FE, DTK_N_TOPO> _finite_elements;
};

template <typename DeviceType>
template <typename Scalar>
Kokkos::View<int *, DeviceType>
Interpolation<DeviceType>::apply( Kokkos::View<Scalar **, DeviceType> X,
                                  Kokkos::View<Scalar **, DeviceType> Y )
{
    // Check that the input and the output have the same number of fields
    DTK_REQUIRE( X.extent( 1 ) == Y.extent( 1 ) );
    using ExecutionSpace = typename DeviceType::execution_space;
    ExecutionSpace space;
    unsigned int const n_fields = X.extent( 1 );
    // Allocate a View that will be used as buffer for the MPI communication
    unsigned int n_local_ref_pts = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        n_local_ref_pts += _point_search._reference_points[topo_id].extent( 0 );
    Kokkos::View<Scalar **, DeviceType> Y_buffer( "Y_buffer", n_local_ref_pts,
                                                  n_fields );

    unsigned int offset = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int const n_ref_points =
            _point_search._reference_points[topo_id].extent( 0 );

        if ( n_ref_points != 0 )
        {
            // Perform the interpolation itself
            Kokkos::View<Scalar **, DeviceType> Y_fe(
                "Y_fe_" + std::to_string( topo_id ), n_ref_points, n_fields );
            interpolateDispatch( _finite_elements[topo_id], topo_id, X, Y_fe );

            // Put Y_fe in the right place in the buffer
            Kokkos::parallel_for(
                DTK_MARK_REGION( "fill_buffer" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
                KOKKOS_LAMBDA( int const i ) {
                    for ( unsigned int j = 0; j < n_fields; ++j )
                        Y_buffer( offset + i, j ) = Y_fe( i, j );
                } );
            Kokkos::fence();
            offset += n_ref_points;
        }
    }

    // Communicate the results, i.e, Y and the associated query ids
    Kokkos::View<unsigned int *, DeviceType> query_ids( "query_ids",
                                                        n_local_ref_pts );
    unsigned int n_copied_pts = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int const size = _point_search._query_ids[topo_id].extent( 0 );
        auto topo_query_ids = _point_search._query_ids[topo_id];
        Kokkos::parallel_for( DTK_MARK_REGION( "query_ids" ),
                              Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                              KOKKOS_LAMBDA( int const i ) {
                                  query_ids( i + n_copied_pts ) =
                                      topo_query_ids( i );
                              } );
        Kokkos::fence();

        n_copied_pts += size;
    }
    unsigned int n_imports =
        _point_search._target_to_source_distributor.getTotalReceiveLength();
    Kokkos::View<unsigned int *, DeviceType> imported_query_ids(
        "imported_query_ids", n_imports );
    Kokkos::View<Scalar **, DeviceType> imported_Y( "imported_Y", n_imports,
                                                    n_fields );
    ArborX::Details::DistributedTreeImpl<DeviceType>::sendAcrossNetwork(
        space, _point_search._target_to_source_distributor, query_ids,
        imported_query_ids );
    ArborX::Details::DistributedTreeImpl<DeviceType>::sendAcrossNetwork(
        space, _point_search._target_to_source_distributor, Y_buffer,
        imported_Y );

    Kokkos::View<int *, DeviceType> found_query_ids( "found_query_ids",
                                                     Y.extent( 0 ) );
    Kokkos::deep_copy( found_query_ids, -1 );

    if ( n_imports != 0 )
    {
        // Because of the MPI communications and the sorting by topologies, all
        // the queries have been reordered. So we put them back in the initial
        // order using the query ids.
        ArborX::Details::DistributedTreeImpl<DeviceType>::sortResults(
            space, imported_query_ids, imported_query_ids, imported_Y );

        // We have finally all the values in the right order but before we can
        // finally return the values we need to filter the data one more time.
        // Some points are correctly found on multiple cells, e.g., point on
        // vertices, so we need to get rid of the duplicates.
        Kokkos::View<unsigned int *, DeviceType> mask( "mask", n_imports );
        Kokkos::deep_copy( mask, 1 );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "compute_mask" ),
            Kokkos::RangePolicy<ExecutionSpace>( 1, n_imports ),
            KOKKOS_LAMBDA( int const i ) {
                if ( imported_query_ids( i - 1 ) == imported_query_ids( i ) )
                    mask( i ) = 0;
            } );
        Kokkos::fence();

        Kokkos::View<unsigned int *, DeviceType> query_offset( "query_offset",
                                                               n_imports );
        ArborX::exclusivePrefixSum( space, mask, query_offset );

        Kokkos::parallel_for(
            DTK_MARK_REGION( "fill_Y" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
            KOKKOS_LAMBDA( int const i ) {
                if ( ( i == 0 ) || ( imported_query_ids( i - 1 ) !=
                                     imported_query_ids( i ) ) )
                {
                    unsigned int k = query_offset( i );
                    for ( unsigned int j = 0; j < n_fields; ++j )
                        Y( k, j ) = imported_Y( i, j );
                    found_query_ids( k ) = imported_query_ids( i );
                }
            } );
        Kokkos::fence();
    }

    return found_query_ids;
}

template <typename DeviceType>
template <typename Scalar, typename FEOpType>
void Interpolation<DeviceType>::interpolate(
    Kokkos::View<Coordinate **, DeviceType> ref_points,
    Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids,
    Kokkos::View<Scalar **, DeviceType> X,
    Kokkos::View<Scalar **, DeviceType> Y_fe )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Functor::Interpolation<Scalar, FEOpType, DeviceType> interpolation_functor(
        _point_search._dim, ref_points, cell_dofs_ids, X, Y_fe );
    Kokkos::parallel_for(
        DTK_MARK_REGION( "interpolate" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, ref_points.extent( 0 ) ),
        interpolation_functor );
}

template <typename DeviceType>
template <typename Scalar, typename FEOpType>
void Interpolation<DeviceType>::hgradInterpolate(
    Kokkos::View<Coordinate **, DeviceType> ref_points,
    Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids,
    Kokkos::View<Scalar **, DeviceType> X,
    Kokkos::View<Scalar **, DeviceType> Y_fe )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Functor::HgradInterpolation<Scalar, FEOpType, DeviceType>
        interpolation_functor( ref_points, cell_dofs_ids, X, Y_fe );
    Kokkos::parallel_for(
        DTK_MARK_REGION( "interpolate" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, ref_points.extent( 0 ) ),
        interpolation_functor );
}

template <typename DeviceType>
template <typename Scalar>
void Interpolation<DeviceType>::interpolateDispatch(
    FE fe, unsigned int topo_id, Kokkos::View<Scalar **, DeviceType> X,
    Kokkos::View<Scalar **, DeviceType> Y_fe )
{

    switch ( fe )
    {
    case FE::HEX_HCURL_1:
    {
        interpolate<Scalar, HEX_HCURL_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::HEX_HDIV_1:
    {
        interpolate<Scalar, HEX_HDIV_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::HEX_HGRAD_1:
    {
        hgradInterpolate<Scalar, HEX_HGRAD_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::HEX_HGRAD_2:
    {
        hgradInterpolate<Scalar, HEX_HGRAD_2::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::PYR_HGRAD_1:
    {
        hgradInterpolate<Scalar, PYR_HGRAD_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::QUAD_HCURL_1:
    {
        interpolate<Scalar, QUAD_HCURL_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::QUAD_HDIV_1:
    {
        interpolate<Scalar, QUAD_HDIV_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::QUAD_HGRAD_1:
    {
        hgradInterpolate<Scalar, QUAD_HGRAD_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::QUAD_HGRAD_2:
    {
        hgradInterpolate<Scalar, QUAD_HGRAD_2::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::TET_HCURL_1:
    {
        interpolate<Scalar, TET_HCURL_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::TET_HDIV_1:
    {
        interpolate<Scalar, TET_HDIV_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::TET_HGRAD_1:
    {
        hgradInterpolate<Scalar, TET_HGRAD_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::TET_HGRAD_2:
    {
        hgradInterpolate<Scalar, TET_HGRAD_2::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::TRI_HGRAD_1:
    {
        hgradInterpolate<Scalar, TRI_HGRAD_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::TRI_HGRAD_2:
    {
        hgradInterpolate<Scalar, TRI_HGRAD_2::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::WEDGE_HGRAD_1:
    {
        hgradInterpolate<Scalar, WEDGE_HGRAD_1::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    case FE::WEDGE_HGRAD_2:
    {
        hgradInterpolate<Scalar, WEDGE_HGRAD_2::feop_type>(
            _point_search._reference_points[topo_id], _dofs_ids[topo_id], X,
            Y_fe );

        break;
    }
    default:
        throw DataTransferKitNotImplementedException();
    }
    Kokkos::fence();
}
} // namespace DataTransferKit

#endif
