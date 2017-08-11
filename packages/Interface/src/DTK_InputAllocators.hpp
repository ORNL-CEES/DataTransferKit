/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
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
    allocateNodeList( const unsigned space_dim, const size_t local_num_nodes );

    // Allocate a bounding volume list.
    static BoundingVolumeList<ViewProperties...>
    allocateBoundingVolumeList( const unsigned space_dim,
                                const size_t local_num_volumes );

    // Allocate a polyhedron list.
    static PolyhedronList<ViewProperties...> allocatePolyhedronList(
        const unsigned space_dim, const size_t local_num_nodes,
        const size_t local_num_faces, const size_t total_face_nodes,
        const size_t local_num_cells, const size_t total_cell_faces );

    // Allocate a cell list.
    static CellList<ViewProperties...>
    allocateCellList( const unsigned space_dim, const size_t local_num_nodes,
                      const size_t local_num_cells,
                      const size_t total_cell_nodes );

    // Allocate a boundary.
    template <class ListType>
    static void allocateBoundary( const size_t local_num_faces,
                                  ListType &list );

    // Allocate an adjacency list.
    template <class ListType>
    static void allocateAdjacencyList( const size_t total_adjacencies,
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

  private:
    // Get the number of cells in a polyhedron list.
    static size_t listNumCells( const PolyhedronList<ViewProperties...> &list );

    // Get the number of cells in a cell list.
    static size_t listNumCells( const CellList<ViewProperties...> &list );
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
