//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_BasicGeometryManager.hpp
 * \author Stuart R. Slattery
 * \brief High-level manager for basic geometries.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BASICGEOMETRYMANAGER_HPP
#define DTK_BASICGEOMETRYMANAGER_HPP

#include <string>

#include "DTK_Types.hpp"
#include "DTK_FunctionSpace.hpp"
#include "DTK_EntitySelector.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class BasicGeometryManager
  \brief High-level manager for basic geometries.

  This manager provides a high-level class for automated construction of DTK
  interface objects. A user is not required to use this class but rather could
  use it to reduce code for certain implementations.
*/
//---------------------------------------------------------------------------//
class BasicGeometryManager
{

    /*!
     * \brief Default constructor. Initializes an empty entity set with a
     * function space defined over the given type.
     */
    BasicGeometryManager( const Teuchos::RCP<const Teuchos::Comm<int> > comm,
			  const int physical_dimension,
			  const EntityType entity_type );

    /*!
     * \brief Entity constructor. Initializes an entity set filled with the
     * specified entities with a function space defined over the given type.
     */
    BasicGeometryManager( const Teuchos::RCP<const Teuchos::Comm<int> > comm,
			  const int physical_dimension,
			  const EntityType entity_type,
			  const Teuchos::ArrayView<Entity>& entities );

    /*!
     * \brief Predicate constructor. Initializes an entity set filled with the
     * specified entities with a selection predicate defined by the input
     * boundary and block ids with a function space defined over the given
     * type.
     */
    BasicGeometryManager( const Teuchos::RCP<const Teuchos::Comm<int> > comm,
			  const int physical_dimension,
			  const EntityType entity_type,
			  const Teuchos::ArrayView<Entity>& entities,
			  const Teuchos::ArrayView<int>& block_ids,
			  const Teuchos::ArrayView<int>& boundary_ids );

    /*!
     * \brief Destructor.
     */
    ~BasicGeometryManager();

    /*!
     * \brief Get the function space over which the mesh and its fields are
     * defined. 
     */
    Teuchos::RCP<FunctionSpace> functionSpace() const;

  private:

    // Create the function space.
    void createFunctionSpace(
	const Teuchos::RCP<const Teuchos::Comm<int> > comm,
	const int physical_dimension, 
	const Teuchos::ArrayView<Entity>& entities,
	const Teuchos::RCP<EntitySelector>& entity_selector );

  private:

    // The function space over which the mesh and its fields are defined.
    Teuchos::RCP<FunctionSpace> d_function_space;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_BASICGEOMETRYMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_BasicGeometryManager.hpp
//---------------------------------------------------------------------------//
