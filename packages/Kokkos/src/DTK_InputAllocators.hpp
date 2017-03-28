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
 * \file DTK_InputAllocators.hpp
 * \brief Allocators for user input data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INPUTALLOCATORS_HPP
#define DTK_INPUTALLOCATORS_HPP

#include "DTK_BoundingVolumeList.hpp"
#include "DTK_CellList.hpp"
#include "DTK_DOFMap.hpp"
#include "DTK_EvaluationSet.hpp"
#include "DTK_Field.hpp"
#include "DTK_NodeList.hpp"
#include "DTK_PolyhedronList.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class InputAllocators
 *
 * \brief Functions for allocating data structures for user input.
 */
template <class... ViewProperties>
class InputAllocators
{
  public:
    // Allocate a node list.
    static NodeList<ViewProperties...>
    allocateNodeList( const unsigned space_dim, const size_t local_num_nodes,
                      const bool has_ghosts );

    // Allocate a bounding volume list.
    static BoundingVolumeList<ViewProperties...>
    allocateBoundingVolumeList( const unsigned space_dim,
                                const size_t local_num_volumes,
                                const bool has_ghosts );

    // Allocate a polyhedron list.
    static PolyhedronList<ViewProperties...> allocatePolyhedronList(
        const unsigned space_dim, const size_t local_num_nodes,
        const size_t local_num_faces, const size_t total_nodes_per_face,
        const size_t local_num_cells, const size_t total_faces_per_cell,
        const bool has_ghosts );

    // Allocate a cell list of cells with the same topology.
    static CellList<ViewProperties...>
    allocateCellList( const unsigned space_dim, const size_t local_num_nodes,
                      const size_t local_num_cells,
                      const unsigned nodes_per_cell, const bool has_ghosts );

    // Allocate a cell list from cells with different topologies.
    static CellList<ViewProperties...> allocateMixedTopologyCellList(
        const unsigned space_dim, const size_t local_num_nodes,
        const size_t local_num_cells, const size_t total_nodes_per_cell,
        const bool has_ghosts );

    // Allocate a boundary.
    template <class ListType>
    static void allocateBoundary( const size_t local_num_faces,
                                  ListType &list );

    // Allocate a degree-of-freedom id Map for objects that all have the same
    // number of degrees of freedom.
    static DOFMap<ViewProperties...>
    allocateDOFMap( const size_t local_num_dofs, const size_t local_num_objects,
                    const unsigned dofs_per_object );

    // Allocate a degree-of-freedom id Map for objects that have
    // different numbers of degrees of freedom.
    static DOFMap<ViewProperties...>
    allocateMixedTopologyDOFMap( const size_t local_num_dofs,
                                 const size_t local_num_objects,
                                 const size_t total_dofs_per_object );
    // Allocate a field.
    template <class Scalar>
    static Field<Scalar, ViewProperties...>
    allocateField( const size_t local_num_dofs,
                   const unsigned field_dimension );

    // Allocate an evaluation set.
    static EvaluationSet<ViewProperties...>
    allocateEvaluationSet( const size_t local_num_evals,
                           const unsigned space_dim );
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Inline implementation.
//---------------------------------------------------------------------------//

#include "DTK_InputAllocators_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_INPUTALLOCATORS_HPP

//---------------------------------------------------------------------------//
// end DTK_InputAllocators.hpp
//---------------------------------------------------------------------------//
