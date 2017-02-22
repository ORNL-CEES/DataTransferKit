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
 * \file DTK_UserApplication.hpp
 * \brief Interface to user applications.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_USERAPPLICATION_HPP
#define DTK_USERAPPLICATION_HPP

#include "DTK_BoundingVolumeList.hpp"
#include "DTK_CellList.hpp"
#include "DTK_ConfigDefs.hpp"
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
    //@{
    //! Type aliases.
    using ExecutionSpace =
        typename ParallelTraits<ParallelModel>::ExecutionSpace;
    //@}

    // Constructor.
    UserApplication(
        const std::shared_ptr<UserFunctionRegistry<Scalar>> &user_functions );

    // Get a node list from the application.
    NodeList<Kokkos::LayoutLeft, ExecutionSpace> getNodeList();

    // Get a bounding volume list from the application.
    BoundingVolumeList<Kokkos::LayoutLeft, ExecutionSpace>
    getBoundingVolumeList();

    // Get a polyhedron list from the application.
    PolyhedronList<Kokkos::LayoutLeft, ExecutionSpace> getPolyhedronList();

    // Get a cell list from the application.
    CellList<Kokkos::LayoutLeft, ExecutionSpace>
    getCellList( std::vector<std::string> &cell_topologies );

    // Get a boundary from the application and put it in the given list.
    template <class ListType>
    void getBoundary( const std::string &boundary_name, ListType &list );

    // Get a dof id map from the application.
    DOFMap<Kokkos::LayoutLeft, ExecutionSpace>
    getDOFMap( std::string &discretization_type );

    // Get a field with a given name from the application.
    Field<Scalar, Kokkos::LayoutLeft, ExecutionSpace>
    getField( const std::string &field_name );

    // Pull a field with a given name to the application.
    void pullField( const std::string &field_name,
                    Field<Scalar, Kokkos::LayoutLeft, ExecutionSpace> field );

    // Push a field with a given name to the application.
    void
    pushField( const std::string &field_name,
               const Field<Scalar, Kokkos::LayoutLeft, ExecutionSpace> field );

    // Ask the application to evaluate a field with a given name.
    void evaluateField(
        const std::string &field_name,
        const EvaluationSet<Kokkos::LayoutLeft, ExecutionSpace> eval_set,
        Field<Scalar, Kokkos::LayoutLeft, ExecutionSpace> field );

  private:
    // User function registry for this application.
    std::shared_ptr<UserFunctionRegistry<Scalar>> _user_functions;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Inline implementation.
//---------------------------------------------------------------------------//

#include "DTK_UserApplication_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_USERAPPLICATION_HPP

//---------------------------------------------------------------------------//
// end DTK_UserApplication.hpp
//---------------------------------------------------------------------------//
