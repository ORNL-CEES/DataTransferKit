/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_POINT_SEARCH_DECL_HPP
#define DTK_POINT_SEARCH_DECL_HPP

#include "DTK_ConfigDefs.hpp"
#include <DTK_CellTypes.h>
#include <DTK_DetailsBox.hpp>
#include <DTK_DetailsPoint.hpp>
#include <DTK_DistributedSearchTree.hpp>
#include <DTK_PointInCell.hpp>

#include <Kokkos_Core.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Distributor.hpp>

#include <tuple>

namespace DataTransferKit
{
template <typename DeviceType>
class PointSearch
{
  public:
    /**
     * Constructor.
     * @param comm
     * @param cell_topologies_view
     * @param cells
     * @param nodes_coordinates coordinates of all the nodes in the mesh
     * @param points_coordinates coordinates in the physical frame of the points
     * that we want to evaluate.
     * @param strategy strategy used when a point is found in multiple cells,
     * e.g., the point is on a face. The two possible strategies are \em unique
     * and \em all. If \em unique is chosen, a single cell will be returned even
     * if the point is found multiple times. If \em all is chosen, all the cells
     * where the point is found will be returned.
     */
    PointSearch( Teuchos::RCP<const Teuchos::Comm<int>> comm,
                 Kokkos::View<DTK_CellTopology *, DeviceType> const
                     &cell_topologies_view,
                 Kokkos::View<unsigned int *, DeviceType> const &cells,
                 Kokkos::View<double **, DeviceType> const &nodes_coordinates,
                 Kokkos::View<double **, DeviceType> const &points_coordinates,
                 std::string const &strategy );

    std::tuple<Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
               Kokkos::View<Point *, DeviceType>,
               Kokkos::View<unsigned int *, DeviceType>>
    getSearchResults();

    void
    convertMesh( std::array<unsigned int, DTK_N_TOPO> const &n_cells_per_topo,
                 Kokkos::View<DTK_CellTopology *, DeviceType> const
                     &cell_topologies_view,
                 Kokkos::View<unsigned int *, DeviceType> const &cells,
                 Kokkos::View<double **, DeviceType> const &coordinates );

    void performDistributedSearch(
        Kokkos::View<double **, DeviceType> points_coord,
        Kokkos::View<Point *, DeviceType> &imported_points,
        Kokkos::View<int *, DeviceType> &imported_query_ids,
        Kokkos::View<int *, DeviceType> &imported_cell_indices,
        Kokkos::View<int *, DeviceType> &ranks );

    template <typename T1, typename T2>
    void computeOffset( Kokkos::View<T1 *, DeviceType> predicate, T2 value,
                        Kokkos::View<unsigned int *, DeviceType> offset );

    void filterTopology( Kokkos::View<unsigned int *, DeviceType> topo,
                         unsigned int topo_id,
                         Kokkos::View<int *, DeviceType> cell_indices,
                         Kokkos::View<Point *, DeviceType> points,
                         Kokkos::View<int *, DeviceType> query_ids,
                         Kokkos::View<int *, DeviceType> ranks,
                         Kokkos::View<unsigned int *, DeviceType> map,
                         Kokkos::View<int *, DeviceType> filtered_cell_indices,
                         Kokkos::View<double **, DeviceType> filtered_points,
                         Kokkos::View<int *, DeviceType> filtered_query_ids,
                         Kokkos::View<int *, DeviceType> filtered_ranks );

    void computeNodeOffset(
        unsigned int const n_cells,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view,
        Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
        Kokkos::View<unsigned int *, DeviceType> nodes_per_cell,
        Kokkos::View<unsigned int *, DeviceType> node_offset );

    void filter_all(
        std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO> const
            &filtered_per_topo_point_in_cell,
        std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO> const
            &filtered_per_topo_reference_points,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
            &filtered_per_topo_cell_indices,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
            &filtered_query_ids,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const &ranks,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
            &filtered_ranks );

    // TODO Interpolation needs to be a friend

  private:
    std::array<unsigned int, DTK_N_TOPO> computeTopologies(
        Kokkos::View<DTK_CellTopology *, DeviceType> const &cell_topo_view );

    void performPointInCell(
        Kokkos::View<double ***, DeviceType> cells,
        Kokkos::View<int *, DeviceType> imported_cell_indices,
        Kokkos::View<Point *, DeviceType> imported_points,
        Kokkos::View<int *, DeviceType> imported_query_ids,
        Kokkos::View<int *, DeviceType> imported_ranks,
        Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
        Kokkos::View<double **, DeviceType> filtered_points,
        Kokkos::View<int *, DeviceType> filtered_cell_indices,
        Kokkos::View<int *, DeviceType> filtered_query_ids,
        Kokkos::View<unsigned int *, DeviceType> points_indices_map,
        Kokkos::View<double **, DeviceType> reference_points,
        Kokkos::View<bool *, DeviceType> point_in_cell,
        Kokkos::View<int *, DeviceType> ranks );

    void filter(
        std::string const &strategy,
        std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO> const
            &filtered_per_topo_point_in_cell,
        std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO> const
            &filtered_per_topo_reference_points,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
            &filtered_per_topo_cell_indices,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
            &filtered_per_topo_query_ids,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const &ranks,
        std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
            &filtered_ranks );

    void build_distributor( std::array<Kokkos::View<int *, DeviceType>,
                                       DTK_N_TOPO> const &filtered_ranks );

    Teuchos::RCP<const Teuchos::Comm<int>> _comm;
    Tpetra::Distributor _target_to_source_distributor;
    unsigned int _dim;
    Kokkos::View<Box *, DeviceType> _bounding_boxes;
    Kokkos::View<unsigned int **, DeviceType> _bounding_box_to_cell;
    Kokkos::View<unsigned int *, DeviceType> _cell_offset;
    std::array<Kokkos::View<Coordinate **, DeviceType>, DTK_N_TOPO>
        _reference_points;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> _query_ids;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> _cell_indices;
    std::array<Kokkos::View<double ***, DeviceType>, DTK_N_TOPO> _block_cells;
};

template <typename DeviceType>
template <typename T1, typename T2>
void PointSearch<DeviceType>::computeOffset(
    Kokkos::View<T1 *, DeviceType> predicate, T2 value,
    Kokkos::View<unsigned int *, DeviceType> offset )
{
    // Create a Kokkos::View with ones where predicate matches value and
    // with zeros everywhere else.
    unsigned int const size = predicate.extent( 0 );
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<unsigned int *, DeviceType> mask( "mask", size );
    Kokkos::parallel_for( REGION_NAME( "compute_mask" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                          KOKKOS_LAMBDA( int const i ) {
                              if ( predicate( i ) == value )
                                  mask( i ) = 1;
                              else
                                  mask( i ) = 0;
                          } );
    Kokkos::fence();

    // Compute an offset that is used be fill filtered_points
    Kokkos::parallel_scan(
        REGION_NAME( "compute_offset" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
        KOKKOS_LAMBDA( int i, int &update, bool final_pass ) {
            if ( final_pass )
                offset( i ) = update;
            update += mask( i );
        } );
    Kokkos::fence();
}
}

#endif
