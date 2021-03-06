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
/*!
 * \file DTK_UserApplication_def.hpp
 * \brief Interface to user applications.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_USERAPPLICATION_DEF_HPP
#define DTK_USERAPPLICATION_DEF_HPP

#include "DTK_InputAllocators.hpp"
#include "DTK_View.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
//! Constructor.
template <class Scalar, class ParallelModel>
UserApplication<Scalar, ParallelModel>::UserApplication(
    const std::shared_ptr<UserFunctionRegistry<Scalar>> &user_functions )
    : _user_functions( user_functions )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Get a node list from the application.
template <class Scalar, class ParallelModel>
auto UserApplication<Scalar, ParallelModel>::getNodeList()
    -> NodeList<Kokkos::LayoutLeft, MemorySpace>
{
    // Get the size of the node list.
    unsigned space_dim;
    size_t local_num_nodes;
    callUserFunction( _user_functions->_node_list_size_func, space_dim,
                      local_num_nodes );

    // Allocate the node list.
    auto node_list =
        InputAllocators<Kokkos::LayoutLeft, MemorySpace>::allocateNodeList(
            space_dim, local_num_nodes );

    // Fill the list with user data.
    View<Coordinate> coordinates( node_list.coordinates );
    callUserFunction( _user_functions->_node_list_data_func, coordinates );

    return node_list;
}

//---------------------------------------------------------------------------//
// Get a bounding volume list from the application.
template <class Scalar, class ParallelModel>
auto UserApplication<Scalar, ParallelModel>::getBoundingVolumeList()
    -> BoundingVolumeList<Kokkos::LayoutLeft, MemorySpace>
{
    // Get the size of the bounding volume list.
    unsigned space_dim;
    size_t local_num_volumes;
    callUserFunction( _user_functions->_bv_list_size_func, space_dim,
                      local_num_volumes );

    // Allocate the bounding volume list.
    auto bv_list = InputAllocators<Kokkos::LayoutLeft, MemorySpace>::
        allocateBoundingVolumeList( space_dim, local_num_volumes );

    // Fill the list with user data.
    View<Coordinate> bounding_volumes( bv_list.bounding_volumes );
    callUserFunction( _user_functions->_bv_list_data_func, bounding_volumes );

    return bv_list;
}

//---------------------------------------------------------------------------//
// Get a polyhedron list from the application.
template <class Scalar, class ParallelModel>
auto UserApplication<Scalar, ParallelModel>::getPolyhedronList()
    -> PolyhedronList<Kokkos::LayoutLeft, MemorySpace>
{
    // Get the size of the polyhedron list.
    unsigned space_dim;
    size_t local_num_nodes;
    size_t local_num_faces;
    size_t total_face_nodes;
    size_t local_num_cells;
    size_t total_cell_faces;
    callUserFunction( _user_functions->_poly_list_size_func, space_dim,
                      local_num_nodes, local_num_faces, total_face_nodes,
                      local_num_cells, total_cell_faces );

    // Allocate the polyhedron list.
    auto poly_list = InputAllocators<Kokkos::LayoutLeft, MemorySpace>::
        allocatePolyhedronList( space_dim, local_num_nodes, local_num_faces,
                                total_face_nodes, local_num_cells,
                                total_cell_faces );

    // Fill the list with user data.
    View<Coordinate> coordinates( poly_list.coordinates );
    View<LocalOrdinal> faces( poly_list.faces );
    View<unsigned> nodes_per_face( poly_list.nodes_per_face );
    View<LocalOrdinal> cells( poly_list.cells );
    View<unsigned> faces_per_cell( poly_list.faces_per_cell );
    View<int> face_orientation( poly_list.face_orientation );
    callUserFunction( _user_functions->_poly_list_data_func, coordinates, faces,
                      nodes_per_face, cells, faces_per_cell, face_orientation );

    return poly_list;
}

//---------------------------------------------------------------------------//
// Get a cell list from the application.
template <class Scalar, class ParallelModel>
auto UserApplication<Scalar, ParallelModel>::getCellList()
    -> CellList<Kokkos::LayoutLeft, MemorySpace>
{
    // Get the size of the cell list.
    unsigned space_dim;
    size_t local_num_nodes;
    size_t local_num_cells;
    size_t total_cell_nodes;
    callUserFunction( _user_functions->_cell_list_size_func, space_dim,
                      local_num_nodes, local_num_cells, total_cell_nodes );

    // Allocate the cell list.
    auto cell_list =
        InputAllocators<Kokkos::LayoutLeft, MemorySpace>::allocateCellList(
            space_dim, local_num_nodes, local_num_cells, total_cell_nodes );

    // Fill the list with user data.
    View<Coordinate> coordinates( cell_list.coordinates );
    View<LocalOrdinal> cells( cell_list.cells );
    View<DTK_CellTopology> cell_topologies( cell_list.cell_topologies );
    callUserFunction( _user_functions->_cell_list_data_func, coordinates, cells,
                      cell_topologies );

    return cell_list;
}

//---------------------------------------------------------------------------//
// Get a boundary from the application.
template <class Scalar, class ParallelModel>
template <class ListType>
void UserApplication<Scalar, ParallelModel>::getBoundary( ListType &list )
{
    // Get the size of the boundary.
    size_t local_num_faces;
    callUserFunction( _user_functions->_boundary_size_func, local_num_faces );

    // Allocate the boundary.
    InputAllocators<Kokkos::LayoutLeft, MemorySpace>::allocateBoundary(
        local_num_faces, list );

    // Fill the boundary with user data.
    View<LocalOrdinal> boundary_cells( list.boundary_cells );
    View<unsigned> cell_faces_on_boundary( list.cell_faces_on_boundary );
    callUserFunction( _user_functions->_boundary_data_func, boundary_cells,
                      cell_faces_on_boundary );
}

//---------------------------------------------------------------------------//
// Get an adjacency list from the application.
template <class Scalar, class ParallelModel>
template <class ListType>
void UserApplication<Scalar, ParallelModel>::getAdjacencyList( ListType &list )
{
    // Get the size of the adjacency list.
    size_t total_adjacencies;
    callUserFunction( _user_functions->_adjacency_list_size_func,
                      total_adjacencies );

    // Allocate the adjacency list.
    InputAllocators<Kokkos::LayoutLeft, MemorySpace>::allocateAdjacencyList(
        total_adjacencies, list );

    // Fill the adjacency list with user data.
    View<GlobalOrdinal> cell_global_ids( list.cell_global_ids );
    View<GlobalOrdinal> adjacent_cell_global_ids( list.adjacent_cells );
    View<unsigned> adjacencies_per_cell( list.adjacencies_per_cell );
    callUserFunction( _user_functions->_adjacency_list_data_func,
                      cell_global_ids, adjacent_cell_global_ids,
                      adjacencies_per_cell );
}

//---------------------------------------------------------------------------//
// Get a dof map from the application.
template <class Scalar, class ParallelModel>
auto UserApplication<Scalar, ParallelModel>::getDOFMap(
    std::string &discretization_type )
    -> DOFMap<Kokkos::LayoutLeft, MemorySpace>
{
    // Both types of dof id maps should not be defined.
    DTK_INSIST( !( _user_functions->_dof_map_size_func.first ) !=
                !( _user_functions->_mt_dof_map_size_func.first ) );

    DOFMap<Kokkos::LayoutLeft, MemorySpace> dof_map;

    // Single topology case.
    if ( _user_functions->_dof_map_size_func.first )
    {
        // Get the size of the dof id map.
        size_t local_num_dofs;
        size_t local_num_objects;
        unsigned dofs_per_object;
        callUserFunction( _user_functions->_dof_map_size_func, local_num_dofs,
                          local_num_objects, dofs_per_object );

        // Allocate the map.
        dof_map =
            InputAllocators<Kokkos::LayoutLeft, MemorySpace>::allocateDOFMap(
                local_num_dofs, local_num_objects, dofs_per_object );

        // Fill the map with user data.
        View<GlobalOrdinal> global_dof_ids( dof_map.global_dof_ids );
        View<LocalOrdinal> object_dof_ids( dof_map.object_dof_ids );
        callUserFunction( _user_functions->_dof_map_data_func, global_dof_ids,
                          object_dof_ids, discretization_type );
    }

    // Multiple topology case.
    else
    {
        // Get the size of the dof id map.
        size_t local_num_dofs;
        size_t local_num_objects;
        size_t total_dofs_per_object;
        callUserFunction( _user_functions->_mt_dof_map_size_func,
                          local_num_dofs, local_num_objects,
                          total_dofs_per_object );

        // Allocate the map.
        dof_map = InputAllocators<Kokkos::LayoutLeft, MemorySpace>::
            allocateMixedTopologyDOFMap( local_num_dofs, local_num_objects,
                                         total_dofs_per_object );

        // Fill the map with user data.
        View<GlobalOrdinal> global_dof_ids( dof_map.global_dof_ids );
        View<LocalOrdinal> object_dof_ids( dof_map.object_dof_ids );
        View<unsigned> dofs_per_object( dof_map.dofs_per_object );
        callUserFunction( _user_functions->_mt_dof_map_data_func,
                          global_dof_ids, object_dof_ids, dofs_per_object,
                          discretization_type );
    }

    return dof_map;
}

//---------------------------------------------------------------------------//
// Get a field with a given name from the application.
template <class Scalar, class ParallelModel>
auto UserApplication<Scalar, ParallelModel>::getField(
    const std::string &field_name )
    -> Field<Scalar, Kokkos::LayoutLeft, MemorySpace>
{
    // Get the size of the field.
    unsigned field_dim;
    size_t local_num_dofs;

    callUserFunction( _user_functions->_field_size_func, field_name, field_dim,
                      local_num_dofs );

    // Allocate the field.
    auto field = InputAllocators<Kokkos::LayoutLeft, MemorySpace>::
        template allocateField<Scalar>( local_num_dofs, field_dim );

    return field;
}

//---------------------------------------------------------------------------//
// Pull a field with a given name to the application.
template <class Scalar, class ParallelModel>
void UserApplication<Scalar, ParallelModel>::pullField(
    const std::string &field_name,
    Field<Scalar, Kokkos::LayoutLeft, MemorySpace> field )
{
    // Get the field from the user.
    View<Scalar> field_dofs( field.dofs );
    callUserFunction( _user_functions->_pull_field_func, field_name,
                      field_dofs );
}

//---------------------------------------------------------------------------//
// Push a field with a given name to the application.
template <class Scalar, class ParallelModel>
void UserApplication<Scalar, ParallelModel>::pushField(
    const std::string &field_name,
    const Field<Scalar, Kokkos::LayoutLeft, MemorySpace> field )
{
    // Give the field to the user.
    View<Scalar> field_dofs( field.dofs );
    callUserFunction( _user_functions->_push_field_func, field_name,
                      field_dofs );
}

//---------------------------------------------------------------------------//
// Ask the application to evaluate a field with a given name.
template <class Scalar, class ParallelModel>
void UserApplication<Scalar, ParallelModel>::evaluateField(
    const std::string &field_name,
    const EvaluationSet<Kokkos::LayoutLeft, MemorySpace> eval_set,
    Field<Scalar, Kokkos::LayoutLeft, MemorySpace> field )
{
    // Ask the user to evaluate the field.
    View<Coordinate> evaluation_points( eval_set.evaluation_points );
    View<LocalOrdinal> object_ids( eval_set.object_ids );
    View<Scalar> values( field.dofs );
    callUserFunction( _user_functions->_eval_field_func, field_name,
                      evaluation_points, object_ids, values );
}

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_USERAPPLICATION_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_UserApplication_def.hpp
//---------------------------------------------------------------------------//
