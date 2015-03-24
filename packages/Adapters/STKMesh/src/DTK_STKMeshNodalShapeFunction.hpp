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
 * \brief DTK_STKMeshNodalShapeFunction.hpp
 * \author Stuart R. Slattery
 * \brief Nodal shape function implementation for STK mesh.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHNODALSHAPEFUNCTION
#define DTK_STKMESHNODALSHAPEFUNCTION

#include "DTK_EntityShapeFunction.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <Intrepid_Basis.hpp>
#include <Intrepid_FieldContainer.hpp>

#include <stk_mesh/base/BulkData.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshNodalShapeFunction
  \brief Nodal shape function implementation for STK mesh.

  STKMeshNodalShapeFunction provides a shape function for node-centered
  quantities with shape functions evaluated in an element supported by
  nodes. The node ids serve as the dof ids for these shape functions. A
  corresponding DOF vector indexed via node ids should be produced to match
  this shape function. STKMeshDOFVector provides services to construct these
  vectors.
*/
//---------------------------------------------------------------------------//
class STKMeshNodalShapeFunction : public EntityShapeFunction
{
  public:

    /*!
     * \brief Constructor.
     */
    STKMeshNodalShapeFunction(
	const Teuchos::RCP<stk::mesh::BulkData>& bulk_data );

    /*!
     * \brief Given an entity, get the ids of the degrees of freedom in the
     * vector space supporting its shape function.
     * \param entity Get the degrees of freedom for this entity.
     * \param dof_ids Return the ids of the degrees of freedom in the parallel
     * vector space supporting the entities.
     */
    void entityDOFIds( const Entity& entity,
		       Teuchos::Array<DofId>& dof_ids ) const;

    /*!
     * \brief Given an entity and a reference point, evaluate the shape
     * function of the entity at that point.
     * \param entity Evaluate the shape function of this entity.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param values Entity shape function evaluated at the reference
     * point. 
     */
    void evaluateValue( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<double> & values ) const;

    /*!
     * \brief Given an entity and a reference point, evaluate the gradient of
     * the shape function of the entity at that point.
     * \param entity Evaluate the shape function of this entity.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param gradients Entity shape function gradients evaluated at the reference
     * point. Return these ordered with respect to those return by
     * getDOFIds() such that gradients[N][D] gives the gradient value of the
     * Nth DOF in the Dth spatial dimension.
     */
    void evaluateGradient( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<Teuchos::Array<double> >& gradients ) const;

  private:

    // Given an entity, get the intrepid basis function.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
    getIntrepidBasis( const Entity& entity ) const;

  private:

    // Bulk data for the mesh over which the shape function is defined.
    Teuchos::RCP<stk::mesh::BulkData> d_bulk_data;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHNODALSHAPEFUNCTION

//---------------------------------------------------------------------------//
// end DTK_STKMeshNodalShapeFunction.hpp
//---------------------------------------------------------------------------//
