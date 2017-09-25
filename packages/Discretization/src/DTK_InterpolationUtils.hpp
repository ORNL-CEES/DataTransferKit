/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_INTERPOLATION_UTILS_HPP
#define DTK_INTERPOLATION_UTILS_HPP

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
template <typename DeviceType>
KOKKOS_FUNCTION void computeBlockCellsBoundingBox(
    unsigned int const dim, int const i, unsigned int const n_nodes,
    unsigned int const node_offset, unsigned int const topo_id,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<unsigned int *, DeviceType> offset,
    Kokkos::View<double **, DeviceType> coordinates,
    Kokkos::View<double ***, DeviceType> block_cells,
    Kokkos::View<DataTransferKit::Box *, DeviceType> bounding_boxes,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell )
{
    DataTransferKit::Box bounding_box;
    // If dim == 2, we need to set bounding_box[4] and
    // bounding_box[5].
    if ( dim == 2 )
    {
        bounding_box[4] = 0;
        bounding_box[5] = 1;
    }
    unsigned int const k = offset( i );
    for ( unsigned int node = 0; node < n_nodes; ++node )
    {
        unsigned int const n = node_offset + node;
        for ( unsigned int d = 0; d < dim; ++d )
        {
            // Copy the coordinated in block_cells
            block_cells( k, node, d ) = coordinates( cells( n ), d );
            // Build the bounding box.
            if ( block_cells( k, node, d ) < bounding_box[d * 2] )
                bounding_box[d * 2] = block_cells( k, node, d );
            if ( block_cells( k, node, d ) > bounding_box[d * 2 + 1] )
                bounding_box[d * 2 + 1] = block_cells( k, node, d );
        }
    }
    bounding_boxes( i ) = bounding_box;
    bounding_box_to_cell( i, topo_id ) = k;
}
}

#endif
