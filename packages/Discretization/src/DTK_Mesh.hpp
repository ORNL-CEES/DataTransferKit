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

#ifndef DTK_MESH_HPP
#define DTK_MESH_HPP

#include <DTK_CellTypes.h>

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
template <typename DeviceType>
struct Mesh
{
  public:
    Mesh( Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_,
          Kokkos::View<unsigned int *, DeviceType> cells_,
          Kokkos::View<Coordinate **, DeviceType> nodes_coordinates_ )
        : cell_topologies( cell_topologies_ )
        , cells( cells_ )
        , nodes_coordinates( nodes_coordinates_ )
    {
    }

    /// cell_topologies (n cells)
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies;
    /// Cells vertices associated to each cell (n cells * n vertices per cell)
    Kokkos::View<unsigned int *, DeviceType> cells;
    /// Nodes_coordinates coordinates of all the nodes in the mesh( n
    /// vertices, dim )
    Kokkos::View<Coordinate **, DeviceType> nodes_coordinates;
};
} // namespace DataTransferKit

#endif
