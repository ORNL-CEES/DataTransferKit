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
 * \brief DTK_EntityIntegrationRule.hpp
 * \author Stuart R. Slattery
 * \brief integration rule interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTEGRATIONRULE_HPP
#define DTK_INTEGRATIONRULE_HPP

#include "DTK_Entity.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityIntegrationRule
  \brief integration rule interface.

  EntityIntegrationRule provides numerical quadrature for entities.
*/
//---------------------------------------------------------------------------//
class EntityIntegrationRule
{
  public:
    /*!
     * \brief Constructor.
     */
    EntityIntegrationRule() { /* ... */}

    /*!
     * \brief Destructor.
     */
    virtual ~EntityIntegrationRule() { /* ... */}

    /*!
     * \brief Given an entity and an integration order, get its integration
     * rule.
     *
     * \param entity Get the integration rule for this entity.
     *
     * \param order Get an integration rule of this order.
     *
     * \param reference_points Return the integration points in the reference
     * frame of the entity in this array. If there are N integration points of
     * topological dimension D then this array is of size
     * reference_points[N][D].
     *
     * \param weights Return the weights of the integration points in this
     * array. If there are N integration points this array is of size
     * weights[N].
     */
    virtual void getIntegrationRule(
        const Entity &entity, const int order,
        Teuchos::Array<Teuchos::Array<double>> &reference_points,
        Teuchos::Array<double> &weights ) const = 0;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTEGRATIONRULE_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityIntegrationRule.hpp
//---------------------------------------------------------------------------//
