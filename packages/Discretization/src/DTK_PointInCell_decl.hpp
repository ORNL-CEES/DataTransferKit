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

#ifndef DTK_POINT_IN_CELL_DECL_HPP
#define DTK_POINT_IN_CELL_DECL_HPP

#include "DTK_ConfigDefs.hpp"
#include <DTK_CellTypes.h>
#include <DTK_DBC.hpp>

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
template <typename DeviceType>
class PointInCell
{
  public:
    /**
     * Performs the local search.
     *    @param[in] physical_points The coordinates of the points in the
     * physical space (coarse_output_size, dim)
     *    @param[in] cells Cells owned by the processor (n_cells, n_nodes, dim)
     *    @param[in] coarse_search_output_cells Indices of local cells from the
     * coarse search (coarse_output_size)
     *    @param[in] cell_topo Topology of the cells in \p cells
     *    @param[out] reference_points The coordinates of the points in the
     * reference space (coarse_output_size, dim)
     *    @param[out] point_in_cell Booleans with value true if the point is in
     * the cell and false otherwise (coarse_output_size)
     */
    static void
    search( Kokkos::View<Coordinate **, DeviceType> physical_points,
            Kokkos::View<Coordinate ***, DeviceType> cells,
            Kokkos::View<int *, DeviceType> coarse_search_output_cells,
            DTK_CellTopology cell_topo,
            Kokkos::View<Coordinate **, DeviceType> reference_points,
            Kokkos::View<bool *, DeviceType> point_in_cell );

    /**
     * Same function as above. However, the function is virtual so that the user
     * can provide their own implementation. If the function is not overriden,
     * it throws an exception.
     *    @param[in] physical_points The coordinates of the points in the
     * physical space (coarse_output_size, dim)
     *    @param[in] cells Cells owned by the processor (n_cells, n_nodes, dim)
     *    @param[in] coarse_search_output_cells Indices of local cells from the
     * coarse search (coarse_output_size)
     *    @param[in] cell_topo Topology of the cells in \p cells
     *    @param[out] reference_points The coordinates of the points in the
     * reference space (coarse_output_size, dim)
     *    @param[out] point_in_cell Booleans with value true if the point is in
     * the cell and false otherwise (coarse_output_size)
     */
    virtual void
    search( Kokkos::View<Coordinate **, DeviceType> physical_points,
            Kokkos::View<Coordinate ***, DeviceType> cells,
            Kokkos::View<int *, DeviceType> coarse_search_output_cells,
            std::string cell_topo,
            Kokkos::View<Coordinate **, DeviceType> reference_points,
            Kokkos::View<bool *, DeviceType> point_in_cell )
    {
        (void)physical_points;
        (void)cells;
        (void)coarse_search_output_cells;
        (void)cell_topo;
        (void)reference_points;
        (void)point_in_cell;

        throw DataTransferKitNotImplementedException();
    }

    static double threshold;
};

// Default value for threshold matches the inclusion tolerance in DTK-2.0 which
// is arbitrary and might need adjustement in client code. See
// https://github.com/ORNL-CEES/DataTransferKit/blob/dtk-2.0/packages/Adapters/Libmesh/src/DTK_LibmeshEntityLocalMap.cpp#L58
template <typename DeviceType>
double PointInCell<DeviceType>::threshold = 1e-6;
} // namespace DataTransferKit

#endif
