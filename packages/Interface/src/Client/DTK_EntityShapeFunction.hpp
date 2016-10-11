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
 * \brief DTK_EntityShapeFunction.hpp
 * \author Stuart R. Slattery
 * \brief Shape function interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SHAPEFUNCTION_HPP
#define DTK_SHAPEFUNCTION_HPP

#include "DTK_Entity.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityShapeFunction
  \brief Shape function interface.

  EntityShapeFunction binds support ids to an entity and provides kernels to
  evaluate the function.
*/
//---------------------------------------------------------------------------//
class EntityShapeFunction
{
  public:
    /*!
     * \brief Constructor.
     */
    EntityShapeFunction();

    /*!
     * \brief Destructor.
     */
    virtual ~EntityShapeFunction();

    /*!
     * \brief Given an entity, get the ids of its support locations.
     *
     * \param entity Get the support ids for this entity.
     *
     * \param support_ids Return the ids of the support locations for the
     * given entity in this array.
     */
    virtual void
    entitySupportIds( const Entity &entity,
                      Teuchos::Array<SupportId> &support_ids ) const = 0;

    /*!
     * \brief Given an entity and a reference point, evaluate shape
     * functions of the entity at that point.
     *
     * \param entity Evaluate shape functions of this entity.
     *
     * \param reference_point Evaluate shape functions at this point
     * given in reference coordinates.
     *
     * \param values Entity shape functions evaluated at the reference
     * point. Return these ordered with respect to those returned by
     * entitySupportIds() such that values[N] gives the value of a shape
     * function of the Nth support location of entity.
     */
    virtual void
    evaluateValue( const Entity &entity,
                   const Teuchos::ArrayView<const double> &reference_point,
                   Teuchos::Array<double> &values ) const = 0;

    /*!
     * \brief Given an entity and a reference point, evaluate the gradient of
     * the shape function of the entity at that point. A default
     * implementation is provided using a finite difference scheme.
     *
     * \param entity Evaluate the shape function of this entity.
     *
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     *
     * \param gradients Entity shape function gradients evaluated at the
     * reference point. Return these ordered with respect to those return by
     * entitySupportIds() such that gradients[N][D] gives the gradient value
     * of the Nth support location in the Dth spatial dimension.
     */
    virtual void
    evaluateGradient( const Entity &entity,
                      const Teuchos::ArrayView<const double> &reference_point,
                      Teuchos::Array<Teuchos::Array<double>> &gradients ) const;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SHAPEFUNCTION_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityShapeFunction.hpp
//---------------------------------------------------------------------------//
