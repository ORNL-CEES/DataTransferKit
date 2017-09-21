/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file DTK_UserDataInterface.hpp
 * \brief C++ DTK interface for accessing application data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_USERDATAINTERFACE_HPP
#define DTK_USERDATAINTERFACE_HPP

#include "DTK_Enumerations.h"
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
/*!
 * \brief User-provided point-in-entity function. Determine if each point is
 * located within the associated candidate entity.
 */
using PointInEntityFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        const View<Coordinate> points,
                        const View<LocalOrdinal> candidate_entity_ids,
                        View<bool> point_in_entity )>;

//---------------------------------------------------------------------------//
// Field interface.
//---------------------------------------------------------------------------//
/*!
 * \brief Get the size parameters for a field layout.
 */
using FieldLayoutSizeFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        size_t& local_num_field_ids,
                        size_t& total_entity_field_ids )>;

//---------------------------------------------------------------------------//
/*
 * \brief Get the data for a field layout.
 */
using FieldLayoutDataFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        View<GlobalOrdinal> local_field_ids,
                        View<GlobalOrdinal> entity_global_field_ids,
                        DTK_Discretization& discretization_type )>;

//---------------------------------------------------------------------------//
/*!
 * \brief Pull data from application into a field. The number of field
 * components is specfied at the time the operator apply function is
 * called. This, along with the field layout, dictates the size of the
 * allocated field values view.
 */
template <class Scalar>
using PullFieldDataFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        View<Scalar> field_values )>;

//---------------------------------------------------------------------------//
/*
 * \brief Push data from a field into the application. The number of field
 * components is specfied at the time the operator apply function is
 * called. This, along with the field layout, dictates the size of the
 * allocated field values view.
 */
template <class Scalar>
using PushFieldDataFunction =
    std::function<void( std::shared_ptr<void> user_data,
                        const View<Scalar> field_values )>;

//---------------------------------------------------------------------------//
// Extended Field Interface
//---------------------------------------------------------------------------//
/*!
 * \brief Assign discretization order to cells and define the basis point type
 * for the discretization. If this function is called the base topology of the
 * cell is assumed to be linear and the n^th order basis function is used. If
 * this function is not defined the order of the cells is assumed to match the
 * topology.
 */
using FieldDiscretizationOrderFunc = std::function<void(
    std::shared_ptr<void> user_data,
    View<int> cell_order,
    DTK_BasisPointType& point_type )>;

//---------------------------------------------------------------------------//
/*
 * \brief User-provided field evaluation function. Evaluate a field at a given
 * set of points in a given set of entities.
 */
template <class Scalar>
using EvaluateFieldFunction = std::function<void(
    std::shared_ptr<void> user_data,
    const View<Coordinate> evaluation_points,
    const View<LocalOrdinal> entity_ids, View<Scalar> field_values )>;

//---------------------------------------------------------------------------//

} // end namespace UserDataInterface

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_USERDATAINTERFACE_HPP

//---------------------------------------------------------------------------//
// end DTK_UserDataInterface.hpp
//---------------------------------------------------------------------------//
