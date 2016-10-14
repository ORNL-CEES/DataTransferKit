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
 * \brief DTK_FunctionSpace.cpp
 * \author Stuart R. Slattery
 * \brief Function space.
 */
//---------------------------------------------------------------------------//

#include "DTK_FunctionSpace.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
//! Constructor.
FunctionSpace::FunctionSpace(
    const Teuchos::RCP<EntitySet> &entity_set,
    const Teuchos::RCP<EntityLocalMap> &local_map,
    const Teuchos::RCP<EntityShapeFunction> &shape_function,
    const Teuchos::RCP<EntityIntegrationRule> &integration_rule,
    const PredicateFunction &select_function )
    : d_entity_set( entity_set )
    , d_local_map( local_map )
    , d_shape_function( shape_function )
    , d_integration_rule( integration_rule )
    , d_select_function( select_function )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Get the entity set over which the fields are defined.
Teuchos::RCP<EntitySet> FunctionSpace::entitySet() const
{
    return d_entity_set;
}

//---------------------------------------------------------------------------//
// Get the reference frame for entities supporting the function.
Teuchos::RCP<EntityLocalMap> FunctionSpace::localMap() const
{
    return d_local_map;
}

//---------------------------------------------------------------------------//
// Get the shape function for entities supporting the function.
Teuchos::RCP<EntityShapeFunction> FunctionSpace::shapeFunction() const
{
    return d_shape_function;
}

//---------------------------------------------------------------------------//
// Get the integration rule for entities supporting the function.
Teuchos::RCP<EntityIntegrationRule> FunctionSpace::integrationRule() const
{
    return d_integration_rule;
}

//---------------------------------------------------------------------------//
// Get the selector function.
PredicateFunction FunctionSpace::selectFunction() const
{
    return d_select_function;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_FunctionSpace.cpp
//---------------------------------------------------------------------------//
