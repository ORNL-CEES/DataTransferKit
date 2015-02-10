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
 * \brief DTK_STKMeshFieldVector_impl.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh field DOF vector.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHFIELDVECTOR_IMPL_HPP
#define DTK_STKMESHFIELDVECTOR_IMPL_HPP

#include <vector>

#include "DTK_DBC.hpp"

#include <Teuchos_DefaultMpiComm.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Scalar, class FieldType>
STKMeshFieldVector<Scalar,FieldType>::STKMeshFieldVector(
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
    const Teuchos::Ptr<FieldType>& field,
    const int field_dim )
    : d_bulk_data( bulk_data )
    , d_field( field )
    , d_field_dim( field_dim )
{
    // Get the field restriction.
    const stk::mesh::FieldRestrictionVector field_restrictions = 
	d_field->restrictions();
    
    // Build a composite selector.
    stk::mesh::Selector field_selector;
    for ( stk::mesh::FieldRestriction r : field_restrictions )
    {
	field_selector = field_selector | r.selector();
    }

    // Get the buckets for the field entity rank.
    const stk::mesh::BucketVector& field_buckets = 
	field_selector.get_buckets( d_field->entity_rank() );
    
    // Get the entities over which the field is defined.
    stk::mesh::get_selected_entities( 
	field_selector, field_buckets, d_field_entities );

    // Extract the field entity ids.
    int num_entities = d_field_entities.size();
    Teuchos::Array<std::size_t> entity_ids( num_entities );
    for ( int n = 0; n < num_entities; ++n )
    {
	entity_ids[n] = d_bulk_data->identifier( d_field_entities[n] );
    }

    // Create a Tpetra map.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::rcp( new Teuchos::MpiComm<int>(d_bulk_data->parallel()) );
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	Tpetra::createNonContigMap<int,std::size_t>( entity_ids, comm );

    // Create a Tpetra vector.
    d_vector =
	Tpetra::createMultiVector<Scalar,int,std::size_t>( map, d_field_dim );
}

//---------------------------------------------------------------------------//
// Given a STK field, create a Tpetra vector that maps to the field DOFs and
// pull the dofs from the field.
template<class Scalar,class FieldType>
void STKMeshFieldVector<Scalar,FieldType>::pullDataFromField()
{
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > vector_view =
	d_vector->get2dViewNonConst();
    Scalar* field_data;
    int num_entity = d_field_entities.size();
    for ( int n = 0; n < num_entity; ++n )
    {
	field_data = stk::mesh::field_data( *d_field, d_field_entities[n] );
	for ( int d = 0; d < d_field_dim; ++d )
	{
	    vector_view[d][n] = field_data[d];
	}
    }
}

//---------------------------------------------------------------------------//
// Given a Tpetra vector of DOF data, push the data into a given STK field.
template<class Scalar,class FieldType>
void STKMeshFieldVector<Scalar,FieldType>::pushDataToField()
{
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > vector_view =
	d_vector->get2dView();
    Scalar* field_data;
    int num_entity = d_field_entities.size();
    for ( int n = 0; n < num_entity; ++n )
    {
	field_data = stk::mesh::field_data( *d_field, d_field_entities[n] );
	for ( int d = 0; d < d_field_dim; ++d )
	{
	    field_data[d] = vector_view[d][n];
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHFIELDVECTOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshFieldVector_impl.hpp
//---------------------------------------------------------------------------//
