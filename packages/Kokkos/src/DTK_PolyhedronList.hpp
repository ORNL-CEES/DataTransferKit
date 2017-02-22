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
 * \file DTK_PolyhedronList.hpp
 * \brief Polyhedron list.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POLYHEDRONLIST_HPP
#define DTK_POLYHEDRONLIST_HPP

#include "DTK_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class PolyhedronList.
 *
 * \brief Trivially-copyable polyhedron list.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 *
 * PolyhedronList describes arbitrary linear polyhedra accessible from the
 * calling MPI rank.
*/
template <class... ViewProperties>
class PolyhedronList
{
  public:
    //! View traits.
    using ViewTraits = typename Kokkos::ViewTraits<int, ViewProperties...>;

    //! The coordinates of the nodes that build the polyhedrons that are
    //! provided through the cell list on this MPI rank. This view is rank-2
    //! and should be sized as (number of nodes, spatial dimension).
    Kokkos::View<Coordinate **, ViewProperties...> coordinates;

    //! Connectivity list of faces that build the polyhedrons that are local
    //! to this MPI rank. This view is defined as rank-1 but represents
    //! unstructured rank-2 data. It should be sized as (total sum of the
    //! number of nodes composing each face) or the sum of all elements in the
    //! following view, nodes_per_face, which indicates how many nodes are
    //! assigned to each face and how to index into this view. The input
    //! should be arranged as follows. Consider the \f$n^th\f$ node of face
    //! \f$i\f$ to be \f$f^i_n\f$ which is equal to the local index of the
    //! corresponding node in the coordinates view. Two faces, the first with
    //! 4 nodes and the second with 3 would then be defined via this view as:
    //! \f$(f^1_1, f^1_2, f^1_3, f^1_4, f^2_1, f^2_2, f^2_3 )\f$ with the
    //! nodes_per_face view reading \f$(4, 3)\f$
    Kokkos::View<LocalOrdinal *, ViewProperties...> faces;

    //! Number of nodes composing each face in the face input view. The sum of
    //! all local elements in this view should equal the total size of the
    //! faces view. This view rank-1 and should be of length of the number of
    //! faces in the list.
    Kokkos::View<unsigned *, ViewProperties...> nodes_per_face;

    //! Connectivity list of polyhedrons. This view is defined as rank-1 but
    //! represents unstructured rank-2 data. It should be sized as (total sum
    //! of the number of faces composing each polyhedron) or the sum of all
    //! elements in the view faces_per_cells, which indicates how many faces
    //! are assigned to each cell and how to index into this view. The input
    //! should be arranged as follows. Consider the \f$n^th\f$ face of cell
    //! \f$i\f$ to be \f$c^i_n\f$ which is equal to the local index of the
    //! corresponding face in the faces view. Two cells, the first with 5
    //! faces and the second with 4 would then be defined via this view as:
    //! \f$(c^1_1, c^1_2, c^1_3, c^1_4, c^1_5, c^2_1, c^2_2, c^2_3, c^2_4 )\f$
    //! with the faces_per_cell view reading \f$(5, 4)\f$.
    Kokkos::View<LocalOrdinal *, ViewProperties...> cells;

    //! Number of cells composing each cell in the cell input view. The sum of
    //! all local elements in this view should equal the total size of the
    //! cells view. This view is rank-1 and of length of the number of cells in
    //! the list.
    Kokkos::View<unsigned *, ViewProperties...> faces_per_cell;

    //! Orientation of each face composing a cell indicating an outward or
    //! inward facing normal based on node ordering of the face and use of the
    //! right-hand rule. This view is defined as rank-1 but represents
    //! unstructured rank-2 data. This view is the same size as the cells view
    //! and is indexed in an identical matter. If the face for the given cell
    //! has a node ordering that returns a face normal that points into the
    //! cell via the right hand rule then a -1 should be input. If the node
    //! ordering of the face produces a normal that points out from the cell a
    //! +1 should be input.
    Kokkos::View<int *, ViewProperties...> face_orientation;

    //! View indicating if the given cell is owned by the local process or is a
    //! ghost. This information is not necessary for all algorithms and
    //! therefore this view can optionally be of size 0 or NULL to indicate
    //! that no ghost information is available thereby indicating that all
    //! cells provided are locally-owned by this MPI rank. This view is rank-1
    //! and of length of the number of cell in the list.
    Kokkos::View<bool *, ViewProperties...> is_ghost_cell;

    //! For every face on the boundary give the local id of the cell to which
    //! the face belongs. This view is of rank-1 and of length equal to the
    //! number of faces on the boundary. If the list does not have a boundary
    //! this view will be empty. Dimensions: (face)
    Kokkos::View<LocalOrdinal *, ViewProperties...> boundary_cells;

    //! For every face on the boundary give the local id of the face relative
    //! to its parent cell. This is the local id of the face relative to the
    //! cell, not the local id of the face in the face view of the polyhedron
    //! list. This view is of rank-1 and of length equal to the number of
    //! faces on the boundary. If the list does not have a boundary this view
    //! will be empty. Dimensions: (face)
    Kokkos::View<unsigned *, ViewProperties...> cell_faces_on_boundary;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_POLYHEDRONLIST_HPP

//---------------------------------------------------------------------------//
// end DTK_PolyhedronList.hpp
//---------------------------------------------------------------------------//
