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
 * \brief DTK_EntityCenteredShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief Shape function interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_EntityCenteredShapeFunction.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
EntityCenteredShapeFunction::EntityCenteredShapeFunction(
    const Teuchos::Array<EntityId>& entity_ids )
    : d_dof_ids( entity_ids )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
EntityCenteredShapeFunction::~EntityCenteredShapeFunction()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get an ordered list of DOF ids on this process. This list must correlate to
// all DOF vectors based on this shape function. This effectively informs us
// of the parallel layout of DOF vectors built on this shape function.
void EntityCenteredShapeFunction::allDOFIds(
    Teuchos::Array<std::size_t>& dof_ids ) const
{
    dof_ids = d_dof_ids;
}

//---------------------------------------------------------------------------//
// Given an entity, get the ids of the degrees of freedom in the vector space
// supporting its shape function.
void EntityCenteredShapeFunction::entityDOFIds( 
    const Entity& entity, Teuchos::Array<std::size_t>& dof_ids ) const
{
    // There is one DOF for an entity-centered quantity. We will assign the
    // DOF id to be the same as the entity id.
    dof_ids.assign( 1, Teuchos::as<std::size_t>(entity.id()) );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void EntityCenteredShapeFunction::evaluateValue( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<double>& values ) const
{
    // There is one DOF and therefore the shape function value is 1.0
    values.assign( 1, 1.0 );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void EntityCenteredShapeFunction::evaluateGradient( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<Teuchos::Array<double> >& gradients ) const
{
    // For now there is no gradient in the entity as we have a constant shape
    // function over its entire domain. In the future, we could use adjacency
    // information to construct some approximation to the derivative via a
    // Taylor expansion.
    gradients.assign( 1, Teuchos::Array<double>(1,0.0) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityCenteredShapeFunction.cpp
//---------------------------------------------------------------------------//
