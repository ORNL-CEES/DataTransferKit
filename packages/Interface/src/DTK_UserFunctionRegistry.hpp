/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file DTK_UserFunctionRegistry.hpp
 * \brief Registry for user functions.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_USERFUNCTIONREGISTRY_HPP
#define DTK_USERFUNCTIONREGISTRY_HPP

#include "DTK_DBC.hpp"
#include "DTK_UserDataInterface.hpp"
#include "DTK_View.hpp"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
//! Use this namespace for convenience.
using namespace UserDataInterface;

//---------------------------------------------------------------------------//
// Forward declaration of UserApplication.
template <class Scalar, class ParallelModel>
class UserApplication;

//---------------------------------------------------------------------------//
/*!
 * \class UserFunctionRegistry
 *
 * \brief Registry for user functions.
 *
 * \tparam Scalar Scalar type of the fields this user function implementation
 * will provide.
 *
 * The registry is the mechanism by which user functions and the data to be
 * called with those user functions are registered with DTK.
 */
//---------------------------------------------------------------------------//
template <class Scalar>
class UserFunctionRegistry
{
    //! UserApplication is a friend class of the registry to provide direct
    //! access to user functions. Only UserApplications with the equivalent
    //! scalar type will access the data of this class. Partial specialization
    //! of friends is prohibited so we fully declare the template parameters
    //! here. Because user functions in the registry are private data and have
    //! no accessors this indicates that the UserApplication class is the only
    //! object allowed to call them.
    template <class UserScalarType, class ParallelModel>
    friend class UserApplication;

  public:
    //! User implementation.
    template <class CallableObject>
    using UserImpl = std::pair<CallableObject, std::shared_ptr<void>>;

    //! @name Set Geometry
    //@{

    //! Node list size function.
    void setNodeListSizeFunction( NodeListSizeFunction &&func,
                                  std::shared_ptr<void> user_data = nullptr );

    //! Node list data function.
    void setNodeListDataFunction( NodeListDataFunction &&func,
                                  std::shared_ptr<void> user_data = nullptr );

    //! Bounding volume list size function.
    void setBoundingVolumeListSizeFunction(
        BoundingVolumeListSizeFunction &&func,
        std::shared_ptr<void> user_data = nullptr );

    //! Bounding volume list data function.
    void setBoundingVolumeListDataFunction(
        BoundingVolumeListDataFunction &&func,
        std::shared_ptr<void> user_data = nullptr );

    //! Polyhedron list size function.
    void
    setPolyhedronListSizeFunction( PolyhedronListSizeFunction &&func,
                                   std::shared_ptr<void> user_data = nullptr );

    //! Polyhedron list data function.
    void
    setPolyhedronListDataFunction( PolyhedronListDataFunction &&func,
                                   std::shared_ptr<void> user_data = nullptr );

    //! Cell list size function.
    void setCellListSizeFunction( CellListSizeFunction &&func,
                                  std::shared_ptr<void> user_data = nullptr );

    //! Cell list data function.
    void setCellListDataFunction( CellListDataFunction &&func,
                                  std::shared_ptr<void> user_data = nullptr );

    //! Mixed topology cell list size function.
    void setMixedTopologyCellListSizeFunction(
        MixedTopologyCellListSizeFunction &&func,
        std::shared_ptr<void> user_data = nullptr );

    //! Mixed topology cell list data function.
    void setMixedTopologyCellListDataFunction(
        MixedTopologyCellListDataFunction &&func,
        std::shared_ptr<void> user_data = nullptr );

    //! Boundary data function.
    void setBoundarySizeFunction( BoundarySizeFunction &&func,
                                  std::shared_ptr<void> user_data = nullptr );

    //! Boundary data function.
    void setBoundaryDataFunction( BoundaryDataFunction &&func,
                                  std::shared_ptr<void> user_data = nullptr );
    //@}

    //! @name Set Degree-of-freedom Maps
    //@{

    //! Single dofs per object dof map size.
    void setDOFMapSizeFunction( DOFMapSizeFunction &&func,
                                std::shared_ptr<void> user_data = nullptr );

    //! Single dofs per object dof map data.
    void setDOFMapDataFunction( DOFMapDataFunction &&func,
                                std::shared_ptr<void> user_data = nullptr );

    //! Multiple dofs per object dof map size.
    void setMixedTopologyDOFMapSizeFunction(
        MixedTopologyDOFMapSizeFunction &&func,
        std::shared_ptr<void> user_data = nullptr );

    //! Multiple dofs per object dof map data.
    void setMixedTopologyDOFMapDataFunction(
        MixedTopologyDOFMapDataFunction &&func,
        std::shared_ptr<void> user_data = nullptr );
    //@}

    //! @name Set Fields
    //@{

    //! Field size.
    void setFieldSizeFunction( FieldSizeFunction<Scalar> &&func,
                               std::shared_ptr<void> user_data = nullptr );

    //! Pull field.
    void setPullFieldDataFunction( PullFieldDataFunction<Scalar> &&func,
                                   std::shared_ptr<void> user_data = nullptr );

    //! Push field.
    void setPushFieldDataFunction( PushFieldDataFunction<Scalar> &&func,
                                   std::shared_ptr<void> user_data = nullptr );

    //! Evaluate field.
    void setEvaluateFieldFunction( EvaluateFieldFunction<Scalar> &&func,
                                   std::shared_ptr<void> user_data = nullptr );
    //@}

  private:
    //@{
    //! User Geometry functions.

    //! Node list size function.
    UserImpl<NodeListSizeFunction> _node_list_size_func;

    //! Node list data function.
    UserImpl<NodeListDataFunction> _node_list_data_func;

    //! Bounding volume size function.
    UserImpl<BoundingVolumeListSizeFunction> _bv_list_size_func;

    //! Bounding volume list data function.
    UserImpl<BoundingVolumeListDataFunction> _bv_list_data_func;

    //! Polyhedron list size function.
    UserImpl<PolyhedronListSizeFunction> _poly_list_size_func;

    //! Polyhedron list data function.
    UserImpl<PolyhedronListDataFunction> _poly_list_data_func;

    //! Single topology cell list size function.
    UserImpl<CellListSizeFunction> _cell_list_size_func;

    //! Single topology cell list data function.
    UserImpl<CellListDataFunction> _cell_list_data_func;

    //! Mixed topology cell list data function.
    UserImpl<MixedTopologyCellListSizeFunction> _mt_cell_list_size_func;

    //! Mixed topology cell list data function.
    UserImpl<MixedTopologyCellListDataFunction> _mt_cell_list_data_func;

    //! Boundary size function.
    UserImpl<BoundarySizeFunction> _boundary_size_func;

    //! Boundary data function.
    UserImpl<BoundaryDataFunction> _boundary_data_func;
    //@}

    //@{
    //! User DOFMap functions.

    //! Single dofs per object dof map size function.
    UserImpl<DOFMapSizeFunction> _dof_map_size_func;

    //! Single dofs per object dof map data function.
    UserImpl<DOFMapDataFunction> _dof_map_data_func;

    //! Multiple dofs per object dof map size function.
    UserImpl<MixedTopologyDOFMapSizeFunction> _mt_dof_map_size_func;

    //! Multiple dofs per object dof map data function.
    UserImpl<MixedTopologyDOFMapDataFunction> _mt_dof_map_data_func;
    //@}

    //@{
    //! User Field functions.

    //! Field size function.
    UserImpl<FieldSizeFunction<Scalar>> _field_size_func;

    //! Field pull data function.
    UserImpl<PullFieldDataFunction<Scalar>> _pull_field_func;

    //! Field push data function.
    UserImpl<PushFieldDataFunction<Scalar>> _push_field_func;

    //! Field evaluate data function.
    UserImpl<EvaluateFieldFunction<Scalar>> _eval_field_func;
    //@}
};

//---------------------------------------------------------------------------//
// Free functions.
//---------------------------------------------------------------------------//
// Call a user function.
template <class UserImpl, class... Args>
void callUserFunction( UserImpl &user_impl, Args &&... args )
{
    DTK_CHECK_USER_FUNCTION( user_impl.first );
    user_impl.first( user_impl.second, std::forward<Args>( args )... );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_USERFUNCTIONREGISTRY_HPP

//---------------------------------------------------------------------------//
// Inline implementation.
//---------------------------------------------------------------------------//

#include "DTK_UserFunctionRegistry_def.hpp"

//---------------------------------------------------------------------------//
// end DTK_UserFunctionRegistry.hpp
//---------------------------------------------------------------------------//
