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
 * \brief DTK_EntityCenteredField_impl.hpp
 * \author Stuart R. Slattery
 * \brief Entity-centered field.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYCENTEREDFIELD_IMPL_HPP
#define DTK_ENTITYCENTEREDFIELD_IMPL_HPP

#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Scalar>
EntityCenteredField<Scalar>::EntityCenteredField(
    const Teuchos::ArrayView<Entity>& entities,
    const int field_dim,
    const Teuchos::ArrayRCP<Scalar>& dof_data,
    const DataLayout layout )
    : d_field_dim( field_dim )
    , d_data( dof_data )
    , d_layout( layout )
{
    d_lda = entities.size();
    DTK_CHECK( dof_data.size() == d_lda * d_field_dim );
	       
    d_dof_ids.resize( d_lda );
    for ( int n = 0; n < d_lda; ++n )
    {
	d_dof_ids[n] = entities[n].id();
	d_id_map.emplace( d_dof_ids[n], n );
    }
}

//---------------------------------------------------------------------------//
// Get the dimension of the field.
template<class Scalar>
int EntityCenteredField<Scalar>::dimension() const
{
    return d_field_dim;
}

//---------------------------------------------------------------------------//
// Get the locally-owned entity DOF ids of the field.
template<class Scalar>
Teuchos::ArrayView<const DofId>
EntityCenteredField<Scalar>::getLocalEntityDOFIds() const
{
    return d_dof_ids();
}

//---------------------------------------------------------------------------//
// Given a local dof id and a dimension, read data from the application
// field.
template<class Scalar>
Scalar EntityCenteredField<Scalar>::readFieldData( const DofId dof_id,
						   const int dimension ) const
{
    DTK_REQUIRE( d_id_map.count(dof_id) );
    int local_id = d_id_map.find( dof_id )->second;
    return (BLOCKED == d_layout) ?
	d_data[dimension*d_lda + local_id] :
	d_data[local_id*d_field_dim + dimension];
}

//---------------------------------------------------------------------------//
// Given a local dof id, dimension, and field value, write data into the
// application field.
template<class Scalar>
void EntityCenteredField<Scalar>::writeFieldData( const DofId dof_id,
						  const int dimension,
						  const Scalar data )
{
    int local_id = d_id_map.find( dof_id )->second;
    switch( d_layout )
    {
	case BLOCKED:
	    d_data[dimension*d_lda + local_id] = data;
	    break;
	case INTERLEAVED:
	    d_data[local_id*d_field_dim + dimension] = data;
	    break;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYCENTEREDFIELD_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityCenteredField_impl.hpp
//---------------------------------------------------------------------------//
