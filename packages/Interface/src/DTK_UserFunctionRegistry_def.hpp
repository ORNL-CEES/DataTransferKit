/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file DTK_UserFunctionRegistry_def.hpp
 * \brief Registry for user functions.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_USERFUNCTIONREGISTRY_DEF_HPP
#define DTK_USERFUNCTIONREGISTRY_DEF_HPP

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Node list size function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setNodeListSizeFunction(
    NodeListSizeFunction &&func, std::shared_ptr<void> user_data )
{
    _node_list_size_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Node list data function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setNodeListDataFunction(
    NodeListDataFunction &&func, std::shared_ptr<void> user_data )
{
    _node_list_data_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Bounding volume list size function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setBoundingVolumeListSizeFunction(
    BoundingVolumeListSizeFunction &&func, std::shared_ptr<void> user_data )
{
    _bv_list_size_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Bounding volume list data function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setBoundingVolumeListDataFunction(
    BoundingVolumeListDataFunction &&func, std::shared_ptr<void> user_data )
{
    _bv_list_data_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Polyhedron list size function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setPolyhedronListSizeFunction(
    PolyhedronListSizeFunction &&func, std::shared_ptr<void> user_data )
{
    _poly_list_size_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Polyhedron list data function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setPolyhedronListDataFunction(
    PolyhedronListDataFunction &&func, std::shared_ptr<void> user_data )
{
    _poly_list_data_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Cell list size function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setCellListSizeFunction(
    CellListSizeFunction &&func, std::shared_ptr<void> user_data )
{
    _cell_list_size_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Cell list data function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setCellListDataFunction(
    CellListDataFunction &&func, std::shared_ptr<void> user_data )
{
    _cell_list_data_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Mixed topology cell list size function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setMixedTopologyCellListSizeFunction(
    MixedTopologyCellListSizeFunction &&func, std::shared_ptr<void> user_data )
{
    _mt_cell_list_size_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Mixed topology cell list data function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setMixedTopologyCellListDataFunction(
    MixedTopologyCellListDataFunction &&func, std::shared_ptr<void> user_data )
{
    _mt_cell_list_data_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Boundary size function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setBoundarySizeFunction(
    const std::string &boundary_name, BoundarySizeFunction &&func,
    std::shared_ptr<void> user_data )
{
    _boundary_size_funcs.emplace( boundary_name,
                                  std::make_pair( func, user_data ) );
}

//---------------------------------------------------------------------------//
// Boundary data function.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setBoundaryDataFunction(
    const std::string &boundary_name, BoundaryDataFunction &&func,
    std::shared_ptr<void> user_data )
{
    _boundary_data_funcs.emplace( boundary_name,
                                  std::make_pair( func, user_data ) );
}

//---------------------------------------------------------------------------//
// Single dofs per object dof map size.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setDOFMapSizeFunction(
    DOFMapSizeFunction &&func, std::shared_ptr<void> user_data )
{
    _dof_map_size_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Single dofs per object dof map data.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setDOFMapDataFunction(
    DOFMapDataFunction &&func, std::shared_ptr<void> user_data )
{
    _dof_map_data_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Multiple dofs per object dof map size.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setMixedTopologyDOFMapSizeFunction(
    MixedTopologyDOFMapSizeFunction &&func, std::shared_ptr<void> user_data )
{
    _mt_dof_map_size_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Multiple dofs per object dof map data.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setMixedTopologyDOFMapDataFunction(
    MixedTopologyDOFMapDataFunction &&func, std::shared_ptr<void> user_data )
{
    _mt_dof_map_data_func = std::make_pair( func, user_data );
}

//---------------------------------------------------------------------------//
// Field size.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setFieldSizeFunction(
    const std::string &field_name, FieldSizeFunction<Scalar> &&func,
    std::shared_ptr<void> user_data )
{
    _field_size_funcs.emplace( field_name, std::make_pair( func, user_data ) );
}

//---------------------------------------------------------------------------//
// Pull field.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setPullFieldDataFunction(
    const std::string &field_name, PullFieldDataFunction<Scalar> &&func,
    std::shared_ptr<void> user_data )
{
    _pull_field_funcs.emplace( field_name, std::make_pair( func, user_data ) );
}

//---------------------------------------------------------------------------//
// Push field.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setPushFieldDataFunction(
    const std::string &field_name, PushFieldDataFunction<Scalar> &&func,
    std::shared_ptr<void> user_data )
{
    _push_field_funcs.emplace( field_name, std::make_pair( func, user_data ) );
}

//---------------------------------------------------------------------------//
// Evaluate field.
template <class Scalar>
void UserFunctionRegistry<Scalar>::setEvaluateFieldFunction(
    const std::string &field_name, EvaluateFieldFunction<Scalar> &&func,
    std::shared_ptr<void> user_data )
{
    _eval_field_funcs.emplace( field_name, std::make_pair( func, user_data ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_USERFUNCTIONREGISTRY_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_UserFunctionRegistry.hpp
//---------------------------------------------------------------------------//
