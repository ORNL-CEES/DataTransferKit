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

#ifndef DTK_POINT_SEARCH_DECL_HPP
#define DTK_POINT_SEARCH_DECL_HPP

#include "DTK_ConfigDefs.hpp"
#include <DTK_CellTypes.h>
#include <DTK_DetailsBox.hpp>
#include <DTK_DetailsPoint.hpp>

#include <Kokkos_Core.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Distributor.hpp>

#include <tuple>

namespace DataTransferKit
{
/**
 * This class performs the search of a set of given points in a given mesh and
 * returns the cell(s) on which each point has been found as well as the
 * position of the points in the reference frame.
 */
template <typename DeviceType>
class PointSearch
{
  public:
    /**
     * Constructor. The search of the points is done in the constructor but
     * the results is not send back to the calling processor.
     * @param comm
     * @param cell_topologies
     * @param cells vertices associated to each cell
     * @param cell_nodes_coordinates coordinates of all the nodes in the mesh
     * @param points_coordinates coordinates in the physical frame of the points
     * that we are looking for.
     * For a more detailed documentation on \p cell_topologies, \p
     * cells, and \p nodes_coordinates see the documentation of CellList.
     */
    PointSearch( Teuchos::RCP<const Teuchos::Comm<int>> comm,
                 Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
                 Kokkos::View<unsigned int *, DeviceType> cells,
                 Kokkos::View<double **, DeviceType> cell_nodes_coordinates,
                 Kokkos::View<double **, DeviceType> points_coordinates );

    /**
     * Return the result of the search. The tuple contains the rank where the
     * points are found, the cell indices associated to the points (local IDs),
     * the coordinates of the points in the frame of reference, and the query
     * ids associated to each point.
     */
    // Note that this function cannot be const because
    // Tpetra::Distributor::doPostsAndWaits is not const
    std::tuple<Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
               Kokkos::View<Point *, DeviceType>,
               Kokkos::View<unsigned int *, DeviceType>>
    getSearchResults();

    /**
     * Create the cells in the format used by Intrepid2.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void buildBlockCells(
        unsigned int n_cells, unsigned int topo_id,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
        Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
        Kokkos::View<unsigned int *, DeviceType> node_offset,
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<unsigned int *, DeviceType> offset,
        Kokkos::View<double **, DeviceType> coordinates );

    /**
     * Build the bounding boxes associated to the cell
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void buildBoundingBoxes(
        unsigned int n_cells, unsigned int topo_id,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
        Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
        Kokkos::View<unsigned int *, DeviceType> node_offset,
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<unsigned int *, DeviceType> offset,
        Kokkos::View<double **, DeviceType> coordinates,
        Kokkos::View<Box *, DeviceType> bounding_boxes );

    /**
     * Build the map between the bounding boxes and the flat array of cells.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void buildBoundingBoxesToBlockCells(
        unsigned int n_cells, unsigned int topo_id,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
        Kokkos::View<unsigned int *, DeviceType> offset,
        Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell );

    /**
     * Convert the 1D Kokkos View cells and coordinates to arrays of 3D Kokkos
     * Views more suitable for Intrepid2.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void convertMesh(
        std::array<unsigned int, DTK_N_TOPO> const &n_cells_per_topo,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<double **, DeviceType> coordinates,
        Kokkos::View<Box *, DeviceType> bounding_boxes,
        Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell );

    /**
     * Perform the distributed search and sends the points and the cell indices
     * to the processors owning the cells.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void performDistributedSearch(
        Kokkos::View<double **, DeviceType> points_coord,
        Kokkos::View<Box *, DeviceType> bounding_boxes,
        Kokkos::View<Point *, DeviceType> &imported_points,
        Kokkos::View<int *, DeviceType> &imported_query_ids,
        Kokkos::View<int *, DeviceType> &imported_cell_indices,
        Kokkos::View<int *, DeviceType> &ranks );

    /**
     * Keep cell_indices, points, query_ids, and ranks that satisfy a given
     * topology.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void filterTopology(
        Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
        Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
        Kokkos::View<int *, DeviceType> cell_indices,
        Kokkos::View<Point *, DeviceType> points,
        Kokkos::View<int *, DeviceType> query_ids,
        Kokkos::View<int *, DeviceType> ranks,
        Kokkos::View<int *, DeviceType> filtered_cell_indices,
        Kokkos::View<double **, DeviceType> filtered_points,
        Kokkos::View<int *, DeviceType> filtered_query_ids,
        Kokkos::View<int *, DeviceType> filtered_ranks );

    /**
     * Keep data corresponding to points found inside the reference cell.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void filterInCell(
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

  private:
    /**
     * Compute the number of cells associated to each topology.
     */
    std::array<unsigned int, DTK_N_TOPO> computeNCellsPerTopology(
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies );

    /**
     * Compute the position in the reference frame of candidates found by the
     * search.
     */
    void performPointInCell(
        Kokkos::View<double ***, DeviceType> cells,
        Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
        Kokkos::View<int *, DeviceType> imported_cell_indices,
        Kokkos::View<Point *, DeviceType> imported_points,
        Kokkos::View<int *, DeviceType> imported_query_ids,
        Kokkos::View<int *, DeviceType> imported_ranks,
        Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
        Kokkos::View<double **, DeviceType> filtered_points,
        Kokkos::View<int *, DeviceType> filtered_cell_indices,
        Kokkos::View<int *, DeviceType> filtered_query_ids,
        Kokkos::View<double **, DeviceType> reference_points,
        Kokkos::View<bool *, DeviceType> point_in_cell,
        Kokkos::View<int *, DeviceType> ranks );

    /**
     * Build the target-to-source distributor.
     */
    void build_distributor( std::array<Kokkos::View<int *, DeviceType>,
                                       DTK_N_TOPO> const &filtered_ranks );

    Teuchos::RCP<const Teuchos::Comm<int>> _comm;
    Tpetra::Distributor _target_to_source_distributor;
    unsigned int _dim;
    std::array<Kokkos::View<Coordinate **, DeviceType>, DTK_N_TOPO>
        _reference_points;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> _query_ids;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> _cell_indices;
    std::array<Kokkos::View<double ***, DeviceType>, DTK_N_TOPO> _block_cells;
    std::array<std::vector<unsigned int>, DTK_N_TOPO> _cell_indices_map;
};
}

#endif
