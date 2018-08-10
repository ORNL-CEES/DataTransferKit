/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file DTK_UserApplication.hpp
 * \brief Interface to user applications.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_USERAPPLICATION_HPP
#define DTK_USERAPPLICATION_HPP

#include "DTK_BoundingVolumeList.hpp"
#include "DTK_CellList.hpp"
#include "DTK_DOFMap.hpp"
#include "DTK_EvaluationSet.hpp"
#include "DTK_Field.hpp"
#include "DTK_NodeList.hpp"
#include "DTK_ParallelTraits.hpp"
#include "DTK_PolyhedronList.hpp"
#include "DTK_UserFunctionRegistry.hpp"
#include "DTK_View.hpp"

#include <Kokkos_Core.hpp>

#include <memory>
#include <string>
#include <vector>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class UserApplication
 *
 * \brief High-level interface to user applications.
 *
 * \tparam Scalar Scalar type of the fields this user application will
 * provide.
 *
 * \tparam ParallelModel The DTK parallel model indicating the on-node
 * parallelism of the user application. Indicates where data will be
 * allocated.
 *
 * The user application provides a high-level interface to compose DTK input
 * data structures and push and pull field data to and from the application
 * through sequences of user function calls.
 */
//---------------------------------------------------------------------------//
template <class Scalar, class ParallelModel>
class UserApplication
{
  public:
    //! @name Type Aliases
    //@{
    using MemorySpace = typename ParallelModel::memory_space;
    //@}

    //! Constructor.
    UserApplication(
        const std::shared_ptr<UserFunctionRegistry<Scalar>> &user_functions );

    //! Get a node list from the application.
    NodeList<Kokkos::LayoutLeft, MemorySpace> getNodeList();

    //! Get a bounding volume list from the application.
    BoundingVolumeList<Kokkos::LayoutLeft, MemorySpace> getBoundingVolumeList();

    //! Get a polyhedron list from the application.
    PolyhedronList<Kokkos::LayoutLeft, MemorySpace> getPolyhedronList();

    //! Get a cell list from the application.
    CellList<Kokkos::LayoutLeft, MemorySpace> getCellList();

    //! Get a boundary from the application and put it in the given list.
    template <class ListType>
    void getBoundary( ListType &list );

    //! Get adjacencies from the application and put it in the given list.
    template <class ListType>
    void getAdjacencyList( ListType &list );

    //! Get a dof id map from the application.
    DOFMap<Kokkos::LayoutLeft, MemorySpace>
    getDOFMap( std::string &discretization_type );

    //! Get a field with a given name from the application.
    Field<Scalar, Kokkos::LayoutLeft, MemorySpace>
    getField( const std::string &field_name );

    //! Pull a field with a given name to the application.
    void pullField( const std::string &field_name,
                    Field<Scalar, Kokkos::LayoutLeft, MemorySpace> field );

    //! Push a field with a given name to the application.
    void
    pushField( const std::string &field_name,
               const Field<Scalar, Kokkos::LayoutLeft, MemorySpace> field );

    //! Ask the application to evaluate a field with a given name.
    void evaluateField(
        const std::string &field_name,
        const EvaluationSet<Kokkos::LayoutLeft, MemorySpace> eval_set,
        Field<Scalar, Kokkos::LayoutLeft, MemorySpace> field );

  private:
    // User function registry for this application.
    std::shared_ptr<UserFunctionRegistry<Scalar>> _user_functions;
};

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//
// Inline implementation.
//---------------------------------------------------------------------------//

#include "DTK_UserApplication_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_USERAPPLICATION_HPP

//---------------------------------------------------------------------------//
// end DTK_UserApplication.hpp
//---------------------------------------------------------------------------//
