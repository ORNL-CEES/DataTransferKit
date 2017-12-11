/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_POINT_IN_CELL_DEF_HPP
#define DTK_POINT_IN_CELL_DEF_HPP

#include <DTK_DBC.hpp>
#include <DTK_PointInCellFunctor.hpp>
#include <DTK_Topology.hpp>
#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
// Because search is static, we cannot use a private function so put the
// function in its own namespace.
namespace internal
{
template <typename CellType, typename DeviceType>
void point_in_cell( double threshold,
                    Kokkos::View<Coordinate **, DeviceType> physical_points,
                    Kokkos::View<Coordinate ***, DeviceType> cells,
                    Kokkos::View<int *, DeviceType> coarse_search_output_cells,
                    Kokkos::View<Coordinate **, DeviceType> reference_points,
                    Kokkos::View<bool *, DeviceType> point_in_cell )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    int const n_ref_pts = reference_points.extent( 0 );

    Functor::PointInCell<CellType, DeviceType> search_functor(
        threshold, physical_points, cells, coarse_search_output_cells,
        reference_points, point_in_cell );
    Kokkos::parallel_for( REGION_NAME( "point_in_cell" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
                          search_functor );
}
}

template <typename DeviceType>
void PointInCell<DeviceType>::search(
    Kokkos::View<Coordinate **, DeviceType> physical_points,
    Kokkos::View<Coordinate ***, DeviceType> cells,
    Kokkos::View<int *, DeviceType> coarse_search_output_cells,
    DTK_CellTopology cell_topo,
    Kokkos::View<Coordinate **, DeviceType> reference_points,
    Kokkos::View<bool *, DeviceType> point_in_cell )
{
    // Check the size of the Views
    DTK_REQUIRE( reference_points.extent( 0 ) == point_in_cell.extent( 0 ) );
    DTK_REQUIRE( reference_points.extent( 0 ) == physical_points.extent( 0 ) );
    DTK_REQUIRE( reference_points.extent( 1 ) == physical_points.extent( 1 ) );
    DTK_REQUIRE( reference_points.extent( 1 ) == cells.extent( 2 ) );

    // Perform the point in cell search. We hide the template parameters used by
    // Intrepid2, using the CellType template.
    // Note that if the Newton solver does not converge, Intrepid2 will just
    // return the last results and there is no way to know that the coordinates
    // in the reference frames where not found.
    if ( cell_topo == DTK_HEX_8 )
    {
        internal::point_in_cell<HEX_8, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_HEX_27 )
    {
        internal::point_in_cell<HEX_27, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_PYRAMID_5 )
    {
        internal::point_in_cell<PYRAMID_5, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_QUAD_4 )
    {
        internal::point_in_cell<QUAD_4, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_QUAD_9 )
    {
        internal::point_in_cell<QUAD_9, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_TET_4 )
    {
        internal::point_in_cell<TET_4, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_TET_10 )
    {
        internal::point_in_cell<TET_10, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_TRI_3 )
    {
        internal::point_in_cell<TRI_3, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_TRI_6 )
    {
        internal::point_in_cell<TRI_6, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_WEDGE_6 )
    {
        internal::point_in_cell<WEDGE_6, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else if ( cell_topo == DTK_WEDGE_18 )
    {
        internal::point_in_cell<WEDGE_18, DeviceType>(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
    }
    else
    {
        throw std::runtime_error( "Not implemented" );
    }
    Kokkos::fence();
}
}

// Explicit instantiation macro
#define DTK_POINTINCELL_INSTANT( NODE )                                        \
    template class PointInCell<typename NODE::device_type>;

#endif
