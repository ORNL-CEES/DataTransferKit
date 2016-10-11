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
 * \brief DTK_FunctionSpace.hpp
 * \author Stuart R. Slattery
 * \brief Function space.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FUNCTIONSPACE_HPP
#define DTK_FUNCTIONSPACE_HPP

#include "DTK_EntityIntegrationRule.hpp"
#include "DTK_EntityLocalMap.hpp"
#include "DTK_EntitySet.hpp"
#include "DTK_EntityShapeFunction.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class FunctionSpace
  \brief Space of a function.

  FunctionSpace binds the functional support of a field to a parallel vector
  space over a selected subset of the entity set.
*/
//---------------------------------------------------------------------------//
class FunctionSpace
{
  public:
    /*!
     * \brief Constructor.
     */
    FunctionSpace( const Teuchos::RCP<EntitySet> &entity_set,
                   const Teuchos::RCP<EntityLocalMap> &local_map,
                   const Teuchos::RCP<EntityShapeFunction> &shape_function,
                   const Teuchos::RCP<EntityIntegrationRule> &integration_rule,
                   const PredicateFunction &select_function = selectAll );

    /*!
     * \brief Get the entity set over which the fields are defined.
     */
    Teuchos::RCP<EntitySet> entitySet() const;

    /*!
     * \brief Get the local map for entities supporting the function.
     */
    Teuchos::RCP<EntityLocalMap> localMap() const;

    /*!
     * \brief Get the shape function for entities supporting the function.
     */
    Teuchos::RCP<EntityShapeFunction> shapeFunction() const;

    /*!
     * \brief Get the integration rule for entities supporting the function.
     */
    Teuchos::RCP<EntityIntegrationRule> integrationRule() const;

    /*!
     * \brief Get the selector function.
     */
    PredicateFunction selectFunction() const;

    /*!
     * \brief Default select function.
     */
    static inline bool selectAll( Entity ) { return true; }

  private:
    // The entity set over which the function space is constructed.
    Teuchos::RCP<EntitySet> d_entity_set;

    // The reference frame for entities in the set.
    Teuchos::RCP<EntityLocalMap> d_local_map;

    // The shape function for the entities in the set.
    Teuchos::RCP<EntityShapeFunction> d_shape_function;

    // The integration rule for the entities in the set.
    Teuchos::RCP<EntityIntegrationRule> d_integration_rule;

    // The selector function.
    PredicateFunction d_select_function;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_FUNCTIONSPACE_HPP

//---------------------------------------------------------------------------//
// end DTK_FunctionSpace.hpp
//---------------------------------------------------------------------------//
