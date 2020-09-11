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
 * \file DTK_UserDataInterface.hpp
 * \brief C++ DTK interface for accessing application data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_USERDATAINTERFACE_HPP
#define DTK_USERDATAINTERFACE_HPP

#include "DTK_CellTypes.h"
#include "DTK_Types.h"
#include "DTK_View.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \namespace UserDataInterface
 *
 * The type alias definitions here serve as argument templates for callable
 * objects to be implemented in the user's application for access to
 * application data by DTK algorithms. The callable objects are meant to be
 * called from the host.
 *
 * The first argument of each callable object is a void shared pointer to user
 * data. This allows the user to register a function with data they would like
 * to use during the implementation of that function. DTK will not access the
 * data in this pointer, rather, it will simply store it when the function is
 * registered and pass it back to the user when the object is called. This is
 * intended to be used in conjunction with std::static_pointer_cast. If the
 * user provided no extra data to be called with the function then this
 * pointer will equal to nullptr when the function is called.
 *
 * If any errors occur during the execution of the user implementation of
 * these functions, a C++ exception should be thrown indicating the error.
 */
//---------------------------------------------------------------------------//

namespace UserDataInterface
{
//---------------------------------------------------------------------------//
// Basic Geometry Interface
//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for building a node list.
 */
using NodeListSizeFunction =
    std::function<void( std::shared_ptr<void> user_data, unsigned &space_dim,
                        size_t &local_num_nodes )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the data for a node list.
 */
using NodeListDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<Coordinate> coordinates )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for building a bounding volume list.
 */
using BoundingVolumeListSizeFunction =
    std::function<void( std::shared_ptr<void> user_data, unsigned &space_dim,
                        size_t &local_num_volumes )>;

//---------------------------------------------------------------------------//
/*
 * \brief Get the data for a bounding volume list.
 */
using BoundingVolumeListDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<Coordinate> bounding_volumes )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for building a polyhedron list.
 */
using PolyhedronListSizeFunction = std::function<void(
    std::shared_ptr<void> user_data, unsigned &space_dim,
    size_t &local_num_nodes, size_t &local_num_faces, size_t &total_face_nodes,
    size_t &local_num_cells, size_t &total_cell_faces )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the data for a polyhedron list.
 */
using PolyhedronListDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<Coordinate> coordinates,
    View<LocalOrdinal> faces, View<unsigned> nodes_per_face,
    View<LocalOrdinal> cells, View<unsigned> faces_per_cell,
    View<int> face_orientation )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for building a cell list.
 */
using CellListSizeFunction =
    std::function<void( std::shared_ptr<void> user_data, unsigned &space_dim,
                        size_t &local_num_nodes, size_t &local_num_cells,
                        size_t &total_cell_nodes )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the data for a cell list.
 */
using CellListDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<Coordinate> coordinates,
    View<LocalOrdinal> cells, View<DTK_CellTopology> cell_topologies )>;

//---------------------------------------------------------------------------//
// Extended Geometry Interface
//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for a boundary.
 */
using BoundarySizeFunction = std::function<void(
    std::shared_ptr<void> user_data, size_t &local_num_faces )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the data for a boundary.
 */
using BoundaryDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<LocalOrdinal> boundary_cells,
    View<unsigned> cell_faces_on_boundary )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for building a cell adjacency list.
 */
using AdjacencyListSizeFunction = std::function<void(
    std::shared_ptr<void> user_data, size_t &total_adjacencies )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the data for an adjacency list.
 */
using AdjacencyListDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<GlobalOrdinal> global_cell_ids,
    View<GlobalOrdinal> adjacent_global_cell_ids,
    View<unsigned> adjacencies_per_cell )>;

//---------------------------------------------------------------------------//
// Degree-of-freedom interface.
//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for a degree-of-freedom id map with a single
 * number of dofs per object.
 */
using DOFMapSizeFunction =
    std::function<void( std::shared_ptr<void> user_data, size_t &local_num_dofs,
                        size_t &local_num_objects, unsigned &dofs_per_object )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the data for a degree-of-freedom id map with a single number of
 * dofs per object.
 */
using DOFMapDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<GlobalOrdinal> global_dof_ids,
    View<LocalOrdinal> object_dof_ids, std::string &discretization_type )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for a degree-of-freedom id map with each
 * object having a potentially different number of dofs (e.g. mixed topology
 * cell lists or polyhedron lists).
 */
using MixedTopologyDOFMapSizeFunction = std::function<void(
    std::shared_ptr<void> user_data, size_t &local_num_dofs,
    size_t &local_num_objects, size_t &total_dofs_per_object )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the data for a multiple object degree-of-freedom id map
 * (e.g. mixed topology cell lists or polyhedron lists).
 */
using MixedTopologyDOFMapDataFunction = std::function<void(
    std::shared_ptr<void> user_data, View<GlobalOrdinal> global_dof_ids,
    View<LocalOrdinal> object_dof_ids, View<unsigned> dofs_per_object,
    std::string &discretization_type )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for a field. Field must be of size
 * local_num_dofs in the associated dof_id_map.
 */
template <class Scalar>
using FieldSizeFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        const std::string &field_name,
                        unsigned &field_dimension, size_t &local_num_dofs )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Pull data from application into a field.
 */
template <class Scalar>
using PullFieldDataFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        const std::string &field_name,
                        View<Scalar> field_dofs )>;

//---------------------------------------------------------------------------//
/*
 * \brief Push data from a field into the application.
 */
template <class Scalar>
using PushFieldDataFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        const std::string &field_name,
                        const View<Scalar> field_dofs )>;

//---------------------------------------------------------------------------//
/*
 * \brief Evaluate a field at a given set of points in a given set of objects.
 */
template <class Scalar>
using EvaluateFieldFunction = std::function<void(
    std::shared_ptr<void> user_data, const std::string &field_name,
    const View<Coordinate> evaluation_points,
    const View<LocalOrdinal> object_ids, View<Scalar> values )>;

//---------------------------------------------------------------------------//

} // namespace UserDataInterface

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_USERDATAINTERFACE_HPP

//---------------------------------------------------------------------------//
// end DTK_UserDataInterface.hpp
//---------------------------------------------------------------------------//
