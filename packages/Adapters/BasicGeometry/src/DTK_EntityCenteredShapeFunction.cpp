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
// Given an entity, get the ids of the degrees of freedom in the vector space
// supporting its shape function.
void EntityCenteredShapeFunction::entitySupportIds(
    const Entity &entity, Teuchos::Array<SupportId> &support_ids ) const
{
    // There is one Support for an entity-centered quantity. We will assign the
    // Support id to be the same as the entity id.
    support_ids.assign( 1, Teuchos::as<SupportId>( entity.id() ) );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void EntityCenteredShapeFunction::evaluateValue(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<double> &values ) const
{
    // There is one Support and therefore the shape function value is 1.0
    values.assign( 1, 1.0 );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void EntityCenteredShapeFunction::evaluateGradient(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<Teuchos::Array<double>> &gradients ) const
{
    // For now there is no gradient in the entity as we have a constant shape
    // function over its entire domain. In the future, we could use adjacency
    // information to construct some approximation to the derivative via a
    // Taylor expansion.
    gradients.assign( 1, Teuchos::Array<double>( 1, 0.0 ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityCenteredShapeFunction.cpp
//---------------------------------------------------------------------------//
