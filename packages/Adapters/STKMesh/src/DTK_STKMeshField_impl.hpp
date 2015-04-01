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
 * \brief DTK_STKMeshField_impl.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh field SUPPORT vector.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHFIELD_IMPL_HPP
#define DTK_STKMESHFIELD_IMPL_HPP

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
STKMeshField<Scalar,FieldType>::STKMeshField(
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
    
    // Get the local entities over which the field is defined.
    stk::mesh::get_selected_entities( 
	field_selector, field_buckets, d_field_entities );

    // Extract the field entity ids.
    int num_entities = d_field_entities.size();
    d_support_ids.resize( num_entities );
    for ( int n = 0; n < num_entities; ++n )
    {
	d_support_ids[n] = d_bulk_data->identifier( d_field_entities[n] );
	d_id_map.emplace( d_support_ids[n], n );
    }
}

//---------------------------------------------------------------------------//
// Get the dimension of the field.
template<class Scalar, class FieldType>
int STKMeshField<Scalar,FieldType>::dimension() const
{
    return d_field_dim;
}

//---------------------------------------------------------------------------//
// Get the locally-owned entity support location ids of the field.
template<class Scalar, class FieldType>
Teuchos::ArrayView<const SupportId>
STKMeshField<Scalar,FieldType>::getLocalSupportIds() const
{
    return d_support_ids();
}

//---------------------------------------------------------------------------//
// Given a local support id and a dimension, read data from the application
// field.
template<class Scalar, class FieldType>
Scalar STKMeshField<Scalar,FieldType>::readFieldData(
    const SupportId support_id, const int dimension ) const
{
    DTK_REQUIRE( d_id_map.count(support_id) );
    int local_id = d_id_map.find( support_id )->second;
    return stk::mesh::field_data(*d_field,d_field_entities[local_id])[dimension];
}

//---------------------------------------------------------------------------//
// Given a local support id, dimension, and field value, write data into the
// application field.
template<class Scalar, class FieldType>
void STKMeshField<Scalar,FieldType>::writeFieldData( const SupportId support_id,
						     const int dimension,
						     const Scalar data )
{
    int local_id = d_id_map.find( support_id )->second;
    stk::mesh::field_data(*d_field,d_field_entities[local_id])[dimension]
	= data;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHFIELD_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshField_impl.hpp
//---------------------------------------------------------------------------//
