/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file DTK_InputAllocators_def.hpp
 * \brief Allocators for user input data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INPUTALLOCATORS_DEF_HPP
#define DTK_INPUTALLOCATORS_DEF_HPP

#include <DTK_CellList.hpp>
#include <DTK_CellTypes.h>
#include <DTK_PolyhedronList.hpp>

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Allocate a node list.
template <class... ViewProperties>
NodeList<ViewProperties...>
InputAllocators<ViewProperties...>::allocateNodeList(
    const unsigned space_dim, const size_t local_num_nodes )
{
    NodeList<ViewProperties...> node_list;

    node_list.coordinates = Kokkos::View<Coordinate **, ViewProperties...>(
        "coordinates", local_num_nodes, space_dim );

    return node_list;
}

//---------------------------------------------------------------------------//
// Allocate a bounding volume list.
template <class... ViewProperties>
BoundingVolumeList<ViewProperties...>
InputAllocators<ViewProperties...>::allocateBoundingVolumeList(
    const unsigned space_dim, const size_t local_num_volumes )
{
    BoundingVolumeList<ViewProperties...> bv_list;

    bv_list.bounding_volumes =
        Kokkos::View<Coordinate * * [2], ViewProperties...>(
            "bounding_volumes", local_num_volumes, space_dim, 2 );

    return bv_list;
}

//---------------------------------------------------------------------------//
// Allocate a polyhedron list.
template <class... ViewProperties>
PolyhedronList<ViewProperties...>
InputAllocators<ViewProperties...>::allocatePolyhedronList(
    const unsigned space_dim, const size_t local_num_nodes,
    const size_t local_num_faces, const size_t total_face_nodes,
    const size_t local_num_cells, const size_t total_cell_faces )
{
    PolyhedronList<ViewProperties...> poly_list;

    poly_list.coordinates = Kokkos::View<Coordinate **, ViewProperties...>(
        "coordinates", local_num_nodes, space_dim );

    poly_list.faces = Kokkos::View<LocalOrdinal *, ViewProperties...>(
        "faces", total_face_nodes );

    poly_list.nodes_per_face = Kokkos::View<unsigned *, ViewProperties...>(
        "nodes_per_face", local_num_faces );

    poly_list.cells = Kokkos::View<LocalOrdinal *, ViewProperties...>(
        "cells", total_cell_faces );

    poly_list.faces_per_cell = Kokkos::View<unsigned *, ViewProperties...>(
        "faces_per_cell", local_num_cells );

    poly_list.face_orientation = Kokkos::View<int *, ViewProperties...>(
        "face_orientation", total_cell_faces );

    return poly_list;
}

//---------------------------------------------------------------------------//
// Allocate a cell list.
template <class... ViewProperties>
CellList<ViewProperties...>
InputAllocators<ViewProperties...>::allocateCellList(
    const unsigned space_dim, const size_t local_num_nodes,
    const size_t local_num_cells, const size_t total_cell_nodes )
{
    CellList<ViewProperties...> cell_list;

    cell_list.coordinates = Kokkos::View<Coordinate **, ViewProperties...>(
        "coordinates", local_num_nodes, space_dim );

    cell_list.cells = Kokkos::View<LocalOrdinal *, ViewProperties...>(
        "cells", total_cell_nodes );

    cell_list.cell_topologies =
        Kokkos::View<DTK_CellTopology *, ViewProperties...>( "cell_topologies",
                                                             local_num_cells );

    return cell_list;
}

//---------------------------------------------------------------------------//
// Allocate a boundary.
template <class... ViewProperties>
template <class ListType>
void InputAllocators<ViewProperties...>::allocateBoundary(
    const size_t local_num_faces, ListType &list )
{
    list.boundary_cells = Kokkos::View<LocalOrdinal *, ViewProperties...>(
        "boundary_cells", local_num_faces );

    list.cell_faces_on_boundary = Kokkos::View<unsigned *, ViewProperties...>(
        "cell_faces_on_boundary", local_num_faces );
}

//---------------------------------------------------------------------------//
// Allocate an adjacency list.
template <class... ViewProperties>
template <class ListType>
void InputAllocators<ViewProperties...>::allocateAdjacencyList(
    const size_t total_adjacencies, ListType &list )
{
    auto num_cells = listNumCells( list );

    list.cell_global_ids = Kokkos::View<GlobalOrdinal *, ViewProperties...>(
        "cell_global_ids", num_cells );

    list.adjacent_cells = Kokkos::View<GlobalOrdinal *, ViewProperties...>(
        "adjacent_cells", total_adjacencies );

    list.adjacencies_per_cell = Kokkos::View<unsigned *, ViewProperties...>(
        "adjacencies_per_cell", num_cells );
}

//---------------------------------------------------------------------------//
// Allocate a degree-of-freedom id Map for objects that all have the same
// number of degrees of freedom.
template <class... ViewProperties>
DOFMap<ViewProperties...> InputAllocators<ViewProperties...>::allocateDOFMap(
    const size_t local_num_dofs, const size_t local_num_objects,
    const unsigned dofs_per_object )
{
    DOFMap<ViewProperties...> dof_id_map;

    dof_id_map.global_dof_ids =
        Kokkos::View<GlobalOrdinal *, ViewProperties...>( "global_dof_ids",
                                                          local_num_dofs );

    dof_id_map.object_dof_ids =
        Kokkos::View<LocalOrdinal **, ViewProperties...>(
            "object_dof_ids", local_num_objects, dofs_per_object );

    return dof_id_map;
}

//---------------------------------------------------------------------------//
// Allocate a degree-of-freedom id Map for objects that have the
// different numbers of degrees of freedom.
template <class... ViewProperties>
DOFMap<ViewProperties...>
InputAllocators<ViewProperties...>::allocateMixedTopologyDOFMap(
    const size_t local_num_dofs, const size_t local_num_objects,
    const size_t total_dofs_per_object )
{
    DOFMap<ViewProperties...> dof_id_map;

    dof_id_map.global_dof_ids =
        Kokkos::View<GlobalOrdinal *, ViewProperties...>( "global_dof_ids",
                                                          local_num_dofs );

    dof_id_map.object_dof_ids = Kokkos::View<LocalOrdinal *, ViewProperties...>(
        "object_dof_ids", total_dofs_per_object );

    dof_id_map.dofs_per_object = Kokkos::View<unsigned *, ViewProperties...>(
        "dofs_per_object", local_num_objects );

    return dof_id_map;
}

//---------------------------------------------------------------------------//
// Allocate a field.
template <class... ViewProperties>
template <class Scalar>
Field<Scalar, ViewProperties...>
InputAllocators<ViewProperties...>::allocateField(
    const size_t local_num_dofs, const unsigned field_dimension )
{
    Field<Scalar, ViewProperties...> field;

    field.dofs = Kokkos::View<Scalar **, ViewProperties...>(
        "dofs", local_num_dofs, field_dimension );

    return field;
}

//---------------------------------------------------------------------------//
// Allocate an evaluation set.
template <class... ViewProperties>
EvaluationSet<ViewProperties...>
InputAllocators<ViewProperties...>::allocateEvaluationSet(
    const size_t local_num_evals, const unsigned space_dim )
{
    EvaluationSet<ViewProperties...> eval_set;

    eval_set.evaluation_points = Kokkos::View<Coordinate **, ViewProperties...>(
        "evaluation_points", local_num_evals, space_dim );

    eval_set.object_ids = Kokkos::View<LocalOrdinal *, ViewProperties...>(
        "object_ids", local_num_evals );

    return eval_set;
}

//---------------------------------------------------------------------------//
// Get the number of cells in a list (PolyhedronList specialization).
template <class... ViewProperties>
size_t InputAllocators<ViewProperties...>::listNumCells(
    const PolyhedronList<ViewProperties...> &list )
{
    return list.faces_per_cell.size();
}

//---------------------------------------------------------------------------//
// Get the number of cells in a list list (CellList specialization).
template <class... ViewProperties>
size_t InputAllocators<ViewProperties...>::listNumCells(
    const CellList<ViewProperties...> &list )
{
    return list.cell_topologies.size();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INPUTALLOCATORS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_InputAllocators.hpp
//---------------------------------------------------------------------------//
