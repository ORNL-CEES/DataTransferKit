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
#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
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

    using ExecutionSpace = typename DeviceType::execution_space;
    int const n_ref_pts = reference_points.extent( 0 );

    // Perform the point in cell search. We hide the template parameters used by
    // Intrepid2, using the CellType template.
    // Note that if the Newton solver does not converge, Intrepid2 will just
    // return the last results and there is no way to know that the coordinates
    // in the reference frames where not found.
    if ( cell_topo == DTK_HEX_8 )
    {
        Functor::PointInCell<CellType::Hexahedron_8, DeviceType> search_functor(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_hex_8" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_HEX_27 )
    {
        Functor::PointInCell<CellType::Hexahedron_27, DeviceType>
            search_functor( threshold, physical_points, cells,
                            coarse_search_output_cells, reference_points,
                            point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_hex_27" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_PYRAMID_5 )
    {
        Functor::PointInCell<CellType::Pyramid_5, DeviceType> search_functor(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_pyr_5" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_QUAD_4 )
    {
        Functor::PointInCell<CellType::Quadrilateral_4, DeviceType>
            search_functor( threshold, physical_points, cells,
                            coarse_search_output_cells, reference_points,
                            point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_quad_4" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_QUAD_9 )
    {
        Functor::PointInCell<CellType::Quadrilateral_9, DeviceType>
            search_functor( threshold, physical_points, cells,
                            coarse_search_output_cells, reference_points,
                            point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_quad_9" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_TET_4 )
    {
        Functor::PointInCell<CellType::Tetrahedron_4, DeviceType>
            search_functor( threshold, physical_points, cells,
                            coarse_search_output_cells, reference_points,
                            point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_tet_4" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_TET_10 )
    {
        Functor::PointInCell<CellType::Tetrahedron_10, DeviceType>
            search_functor( threshold, physical_points, cells,
                            coarse_search_output_cells, reference_points,
                            point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_tet_10" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_TRI_3 )
    {
        Functor::PointInCell<CellType::Triangle_3, DeviceType> search_functor(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_tri_3" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_TRI_6 )
    {
        Functor::PointInCell<CellType::Triangle_6, DeviceType> search_functor(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_tri_6" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_WEDGE_6 )
    {
        Functor::PointInCell<CellType::Wedge_6, DeviceType> search_functor(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_wedge_6" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
    }
    else if ( cell_topo == DTK_WEDGE_18 )
    {
        Functor::PointInCell<CellType::Wedge_18, DeviceType> search_functor(
            threshold, physical_points, cells, coarse_search_output_cells,
            reference_points, point_in_cell );
        Kokkos::parallel_for(
            REGION_NAME( "compute_pos_in_ref_space_wedge_18" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_pts ),
            search_functor );
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
