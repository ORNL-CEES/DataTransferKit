/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file DTK_CellList.hpp
 * \brief Cell list.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CELLLIST_HPP
#define DTK_CELLLIST_HPP

#include "DTK_CellTypes.h"

#include <Kokkos_Core.hpp>

#include <vector>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class CellList
 *
 * \brief Trivially-copyable cell list.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 *
 * CellList describes cells for various discretizations (FEM, FD, FV, etc.)
 * using a standard set of cell topologies indicated via a unique cell
 * topology key.
*/
template <class... ViewProperties>
class CellList
{
  public:
    //! View traits.
    using ViewTraits = typename Kokkos::ViewTraits<int, ViewProperties...>;

    //! The coordinates of the nodes that build the cells that are provided
    //! through the cell list on this MPI rank. This view is rank-2 and should
    //! be sized as (number of nodes, spatial dimension).
    Kokkos::View<Coordinate **, ViewProperties...> coordinates;

    //! Connectivity list of cells. The view rank is rank-1.
    //!
    //! It represents a lists of cells with different topologies ordered in
    //! blocks. It should be sized as (total sum of the number of nodes
    //! composing each cell). The input should be arranged as
    //! follows. Consider the \f$n^th\f$ node of cell \f$i\f$ to be
    //! \f$c^i_n\f$ which is equal to the local index of the corresponding
    //! node in the nodes view. Two cells, the first with 5 nodes and the
    //! second with 4 would then be defined via this view as: \f$(c^1_1,
    //! c^1_2, c^1_3, c^1_4, c^1_5, c^2_1, c^2_2, c^2_3, c^2_4 )\f$ with the
    //! nodes_per_cell view reading \f$(5, 4)\f$. The number of nodes per cell
    //! is defined by the topology of the cell block given by the associated
    //! entry in block_topologies.
    Kokkos::View<LocalOrdinal *, ViewProperties...> cells;

    //! The topologies of the cells.
    Kokkos::View<DTK_CellTopology *, ViewProperties...> cell_topologies;

    //! View indicating if the given cell is owned by the local process or is a
    //! ghost. This information is not necessary for all algorithms and
    //! therefore this view can optionally be of size 0 or NULL to indicate
    //! that no ghost information is available thereby indicating that all
    //! cells provided are locally-owned by this MPI rank. This view is of
    //! rank-1 and of the length of the number of cells in the list.
    Kokkos::View<bool *, ViewProperties...> is_ghost_cell;

    //! For every face on the boundary give the local id of the cell to which
    //! the face belongs. This view is of rank-1 and of length equal to the
    //! number of faces on the boundary. If the list does not have a boundary
    //! this view will be empty. Dimensions: (face)
    Kokkos::View<LocalOrdinal *, ViewProperties...> boundary_cells;

    //! For every face on the boundary give the local id of the face relative
    //! to its parent cell. This is the local face id relative to the nodes as
    //! defined by the canonical cell topology. This view is of rank-1 and of
    //! length equal to the number of faces on the boundary. If the list does
    //! not have a boundary this view will be empty. Dimensions: (face)
    Kokkos::View<unsigned *, ViewProperties...> cell_faces_on_boundary;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CELLLIST_HPP

//---------------------------------------------------------------------------//
// end DTK_CellList.hpp
//---------------------------------------------------------------------------//
