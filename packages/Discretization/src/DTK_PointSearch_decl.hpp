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

#ifndef DTK_POINT_SEARCH_DECL_HPP
#define DTK_POINT_SEARCH_DECL_HPP

#include "DTK_ConfigDefs.hpp"
#include <ArborX.hpp>
#include <DTK_CellTypes.h>
#include <DTK_Mesh.hpp>

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <array>
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
     * @param mesh mesh of the domain of interest
     * @param points_coordinates coordinates in the physical frame of the points
     * that we are looking for.
     * For a more detailed documentation on \p cell_topologies, \p
     * cells, and \p nodes_coordinates see the documentation of CellList.
     */
    PointSearch( MPI_Comm comm, Mesh<DeviceType> const &mesh,
                 Kokkos::View<Coordinate **, DeviceType> points_coordinates );

    /**
     * Return the result of the search. The tuple contains the rank where the
     * points are found, the cell indices associated to the points (local IDs),
     * the coordinates of the points in the frame of reference, and the query
     * ids associated to each point.
     */
    std::tuple<Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
               Kokkos::View<Coordinate * [3], DeviceType>,
               Kokkos::View<unsigned int *, DeviceType>>
    getSearchResults() const;

    /**
     * Perform the distributed search and sends the points and the cell indices
     * to the processors owning the cells.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    std::tuple<Kokkos::View<ArborX::Point *, DeviceType>,
               Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
               Kokkos::View<int *, DeviceType>>
    performDistributedSearch(
        Kokkos::View<Coordinate **, DeviceType> points_coord,
        Kokkos::View<ArborX::Box *, DeviceType> bounding_boxes );

    /**
     * Keep cell_indices, points, query_ids, and ranks that satisfy a given
     * topology.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    std::tuple<Kokkos::View<int *, DeviceType>,
               Kokkos::View<Coordinate **, DeviceType>,
               Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>>
    filterTopology(
        Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
        unsigned int size,
        Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
        Kokkos::View<int *, DeviceType> cell_indices,
        Kokkos::View<ArborX::Point *, DeviceType> points,
        Kokkos::View<int *, DeviceType> query_ids,
        Kokkos::View<int *, DeviceType> ranks );

    /**
     * Keep data corresponding to points found inside the reference cell.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    Kokkos::View<int *, DeviceType> filterInCell(
        Kokkos::View<bool *, DeviceType> filtered_per_topo_point_in_cell,
        Kokkos::View<Coordinate **, DeviceType>
            filtered_per_topo_reference_points,
        Kokkos::View<int *, DeviceType> filtered_per_topo_cell_indices,
        Kokkos::View<int *, DeviceType> filtered_per_topo_query_ids,
        Kokkos::View<int *, DeviceType> filtered_per_topo_ranks,
        unsigned int topo_id );

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
    Kokkos::View<int *, DeviceType> performPointInCell(
        Kokkos::View<Coordinate ***, DeviceType> cells,
        Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
        Kokkos::View<int *, DeviceType> imported_cell_indices,
        Kokkos::View<ArborX::Point *, DeviceType> imported_points,
        Kokkos::View<int *, DeviceType> imported_query_ids,
        Kokkos::View<int *, DeviceType> imported_ranks,
        Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
        unsigned int size );

    /**
     * Build the target-to-source distributor.
     */
    void build_distributor( std::array<Kokkos::View<int *, DeviceType>,
                                       DTK_N_TOPO> const &filtered_ranks );

    template <typename T>
    friend class Interpolation;

    MPI_Comm _comm;
    ArborX::Details::Distributor<DeviceType> _target_to_source_distributor;
    unsigned int _dim;
    std::array<Kokkos::View<Coordinate **, DeviceType>, DTK_N_TOPO>
        _reference_points;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> _query_ids;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> _cell_indices;
    std::array<std::vector<unsigned int>, DTK_N_TOPO> _cell_indices_map;
};
} // namespace DataTransferKit

#endif
