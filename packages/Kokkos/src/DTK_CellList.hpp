//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file DTK_CellList.hpp
 * \brief Cell list.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CELLLIST_HPP
#define DTK_CELLLIST_HPP

#include "DTK_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>

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

    //! Connectivity list of cells. The view rank is defined dynamically and
    //! can be either rank-1 or rank-2.
    //!
    //! If rank-1, it represents an unordered lists of cells with different
    //! topologies. It should be sized as (total sum of the number of nodes
    //! composing each cell). The input should be arranged as follows. Consider
    //! the \f$n^th\f$ node of cell \f$i\f$ to be \f$c^i_n\f$ which is equal to
    //! the local
    //! index of the corresponding node in the nodes view. Two cells, the first
    //! with 5 nodes and the second with 4 would then be defined via this view
    //! as: \f$(c^1_1, c^1_2, c^1_3, c^1_4, c^1_5, c^2_1, c^2_2, c^2_3, c^2_4
    //! )\f$
    //! with the nodes_per_cell view reading \f$(5, 4)\f$. The number of nodes
    //! per
    //! cell is defined by the topology of the cell given by the associated
    //! entry in topology_keys.
    //!
    //! If rank-2, it represents a list of cells all with the same
    //! topologies. It should be dimensioned as (cells, nodes per cell). The
    //! nodes for each cell should be given in the same canonical order as the
    //! rank-1 input case.
    Kokkos::DynRankView<LocalOrdinal, ViewProperties...> cells;

    //! The local topology id for each cell in the cells view. This view is of
    //! rank-1 and of the length of the number of cells in the list. This
    //! should only be filled by the user if the cells view is of rank-1,
    //! indicated potentially a different topology for every cell in the list.
    Kokkos::View<unsigned *, ViewProperties...> cell_topology_ids;

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
