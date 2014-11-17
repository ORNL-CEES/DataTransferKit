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
 * \brief DTK_STKMeshDOFVector_impl.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh field DOF vector.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHDOFVECTOR_IMPL_HPP
#define DTK_STKMESHDOFVECTOR_IMPL_HPP

#include <vector>

#include "DTK_DBC.hpp"

#include <Teuchos_DefaultMpiComm.hpp>

#include <Tpetra_Map.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Given a STK field, create a Tpetra vector that maps to the field DOFs.
template<class Scalar,class FieldType>
Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
STKMeshDOFVector::createTpetraMultiVectorFromSTKField( 
    const stk::mesh::BulkData& bulk_data,
    const FieldType& field,
    const int field_dim )
{
    // Get the field restriction.
    const stk::mesh::FieldRestrictionVector field_restrictions = 
	field.restrictions();
    
    // Build a composite selector.
    stk::mesh::Selector field_selector;
    for ( stk::mesh::FieldRestriction r : field_restrictions )
    {
	field_selector = field_selector | r.selector();
    }

    // Get the buckets for the field entity rank.
    const stk::mesh::BucketVector& field_buckets = 
	field_selector.get_buckets( field.entity_rank() );
    
    // Get the entities over which the field is defined.
    std::vector<stk::mesh::Entity> field_entities;
    stk::mesh::get_selected_entities( 
	field_selector, field_buckets, field_entities );

    // Extract the field entity ids.
    int num_entities = field_entities.size();
    Teuchos::Array<std::size_t> entity_ids( num_entities );
    for ( int n = 0; n < num_entities; ++n )
    {
	entity_ids[n] = bulk_data.identifier( field_entities[n] );
    }

    // Create a Tpetra map.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::rcp( new Teuchos::MpiComm<int>(bulk_data.parallel()) );
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	Tpetra::createNonContigMap<int,std::size_t>( entity_ids, comm );

    // Create a Tpetra vector.
    return Tpetra::createMultiVector<Scalar,int,std::size_t>( map, field_dim );
}

//---------------------------------------------------------------------------//
// Given a Tpetra vector of DOF data, push the data into a given STK field.
template<class Scalar,class FieldType>
void STKMeshDOFVector::pushTpetraMultiVectorToSTKField(
    const Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> >& field_dofs,
    const stk::mesh::BulkData& bulk_data,
    FieldType& field )
{
    // Get a view to the vector data.
    int field_dim = field_dofs->getNumVectors();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > dofs =
		      field_dofs->get2dView();

    // Get the ids of the entities over which the field is defined.
    Teuchos::ArrayView<const std::size_t> field_entity_ids =
	field_dofs->getMap()->getNodeElementList();

    // Extract the field data from the Tpetra vector into the STK mesh.
    Scalar* field_data;
    int num_entity = field_entity_ids.size();
    for ( int n = 0; n < num_entity; ++n )
    {
	const stk::mesh::Entity stk_entity = bulk_data.get_entity( 
	    field.entity_rank(), field_entity_ids[n] );

	field_data = stk::mesh::field_data( field, stk_entity );

	for ( int d = 0; d < field_dim; ++d )
	{
	    field_data[d] = dofs[d][n];
	}
    }
}

//---------------------------------------------------------------------------//
template<class Scalar>
Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
STKMeshDOFVector::createTpetraMultiVectorFromView( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::ArrayView<const std::size_t>& entity_ids,
    const Teuchos::ArrayRCP<Scalar>& dof_data,
    const std::size_t lda,
    const std::size_t num_vectors )
{
    DTK_REQUIRE( lda*num_vectors == Teuchos::as<std::size_t>(dof_data.size()) );

    // Construct a map.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	Tpetra::createNonContigMap<int,std::size_t>( entity_ids, comm );

    // Build a tpetra multivector.
    return Tpetra::createMultiVectorFromView<Scalar,int,std::size_t>( 
	    map, dof_data, lda, num_vectors );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHDOFVECTOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshDOFVector_impl.hpp
//---------------------------------------------------------------------------//
