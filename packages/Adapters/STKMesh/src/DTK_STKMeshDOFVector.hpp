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
 * \brief DTK_STKMeshDOFVector.hpp
 * \author Stuart R. Slattery
 * \brief Helper functions for managing STK mesh DOF vectors.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHDOFVECTOR_HPP
#define DTK_STKMESHDOFVECTOR_HPP

#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_MultiVector.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshDOFVector
  \brief Helper functions for managing STK mesh DOF vectors.
*/
//---------------------------------------------------------------------------//
class STKMeshDOFVector
{
  public:

    /*!
     * \brief Constructor.
     */
    STKMeshDOFVector()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~STKMeshDOFVector()
    { /* ... */ }

    /*!
     * \brief Given a STK field, create a Tpetra vector that maps to the field
     * DOFs.
     * \param bulk_data The bulk data over which the STK field is defined.
     * \param field The STK field.
     * \param field_dim The dimension of the field. This is the product of the
     * size of all ranks.
     * \return A Tpetra MultiVector indexed according to the field entities
     * with a vector for each field dimension.
     */
    template<class Scalar,class FieldType>
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
    createTpetraMultiVectorFromSTKField( const stk::mesh::BulkData& bulk_data,
					 const FieldType& field,
					 const int field_dim );

    /*!
     * \brief Given a Tpetra vector of DOF data, push the data into a given
     * STK field.
     * \param field_dofs A Tpetra vector containing the field DOFs. One vector
     * for each dimension of the field.
     * \param bulk_data The bulk data over which the STK field is defined.
     * \param field The STK field.
     */
    template<class Scalar,class FieldType>
    static void pushTpetraMultiVectorToSTKField(
	const Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> >& field_dofs,
	const stk::mesh::BulkData& bulk_data,
	FieldType& field );

    /*!
     * \brief Given a set of entities and DOF data bound to those entites,
     * build a Tpetra vector.
     * \param bulk_data The bulk data that owns the entities.
     * \param entities The entities the DOFs are defined over.
     * \param field_dim The dimension of the field. This is the product of the
     * size of all ranks.
     * \param dof_data Reference counted array of the DOF data for the given
     * entity ids.
     * \return A Tpetra MultiVector indexed according to the field entities
     * with a vector for each field dimension.
     */
    template<class Scalar>
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
    createTpetraMultiVectorFromEntitiesAndView(
	const stk::mesh::BulkData& bulk_data,
	const std::vector<stk::mesh::Entity>& entities,
	const int field_dim,
	const Teuchos::ArrayRCP<Scalar>& dof_data );

    /*!
     * \brief Given a set of entity ids and DOF data bound to those entites,
     * build a Tpetra vector.
     * \param bulk_data The bulk data that owns the entities.
     * \param parts The parts over which to build the field.
     * \param field_entity_rank The rank of the entities in the part over
     * which to build the field.
     * \param field_dim The dimension of the field. This is the product of the
     * size of all ranks.
     * \param dof_data Reference counted array of the DOF data for the given
     * entity ids.
     * \return A Tpetra MultiVector indexed according to the field entities
     * with a vector for each field dimension.
     */
    template<class Scalar>
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
    createTpetraMultiVectorFromPartVectorAndView(
	const stk::mesh::BulkData& bulk_data,
	const stk::mesh::PartVector& parts,
	const stk::mesh::EntityRank field_entity_rank,
	const int field_dim,
	const Teuchos::ArrayRCP<Scalar>& dof_data );

    /*!
     * \brief Given a set of entity ids and DOF data bound to those entites,
     * build a Tpetra vector.
     * \param bulk_data The bulk data that owns the entities.
     * \param selector The selector over which to build the field.
     * \param field_entity_rank The rank of the entities in the part over
     * which to build the field.
     * \param field_dim The dimension of the field. This is the product of the
     * size of all ranks.
     * \param dof_data Reference counted array of the DOF data for the given
     * entity ids.
     * \return A Tpetra MultiVector indexed according to the field entities
     * with a vector for each field dimension.
     */
    template<class Scalar>
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
    createTpetraMultiVectorFromSelectorAndView(
	const stk::mesh::BulkData& bulk_data,
	const stk::mesh::Selector& selector,
	const stk::mesh::EntityRank field_entity_rank,
	const int field_dim,
	const Teuchos::ArrayRCP<Scalar>& dof_data );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_STKMeshDOFVector_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHDOFVECTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshDOFVector.hpp
//---------------------------------------------------------------------------//
