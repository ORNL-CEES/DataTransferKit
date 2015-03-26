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
 * \brief DTK_MoabTagField_impl.hpp
 * \author Stuart R. Slattery
 * \brief Moab tag vector manager
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABTAGFIELD_IMPL_HPP
#define DTK_MOABTAGFIELD_IMPL_HPP

#include <vector>

#include "DTK_DBC.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Scalar>
MoabTagField<Scalar>::MoabTagField(
    const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
    const moab::EntityHandle& mesh_set,
    const moab::Tag& tag )
    : d_moab_mesh( moab_mesh )
    , d_mesh_set( mesh_set )
    , d_tag( tag )
{
    // Get the dimension of the tag.
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->tag_get_length( d_tag, d_tag_dim )
	);

    // Get the entities in the set.
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->get_entities_by_handle( d_mesh_set, d_entities );
	);

    // Create local ids.
    int num_entities = d_entities.size();
    d_dof_ids.resize( num_entities );
    for ( int n = 0; n < num_entities; ++n )
    {
	d_dof_ids[n] = d_entities[n];
	d_id_map.emplace( d_dof_ids[n], n );
    }
}

//---------------------------------------------------------------------------//
// Get the dimension of the field.
template<class Scalar>
int MoabTagField<Scalar>::dimension() const
{
    return d_tag_dim;
}

//---------------------------------------------------------------------------//
// Get the locally-owned entity DOF ids of the field.
template<class Scalar>
Teuchos::ArrayView<const DofId>
MoabTagField<Scalar>::getLocalEntityDOFIds() const
{
    return d_dof_ids();
}

//---------------------------------------------------------------------------//
// Given a local dof id and a dimension, read data from the application
// field.
template<class Scalar>
Scalar MoabTagField<Scalar>::readFieldData( const DofId dof_id,
					    const int dimension ) const
{
    DTK_REQUIRE( d_id_map.count(dof_id) );
    int local_id = d_id_map.find( dof_id )->second;
    const void* tag_data = 0;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->tag_get_by_ptr(
	    d_tag,
	    &d_entities[local_id],
	    1,
	    &tag_data )
	);
    return static_cast<const Scalar*>(tag_data)[dimension];
}

//---------------------------------------------------------------------------//
// Given a local dof id, dimension, and field value, write data into the
// application field.
template<class Scalar>
void MoabTagField<Scalar>::writeFieldData( const DofId dof_id,
					   const int dimension,
					   const Scalar data )
{
    int local_id = d_id_map.find( dof_id )->second;
    const void* tag_data = 0;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->tag_get_by_ptr(
	    d_tag,
	    &d_entities[local_id],
	    1,
	    &tag_data )
	);
    const_cast<Scalar*>(static_cast<const Scalar*>(tag_data))[dimension] = data;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABTAGFIELD_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabTagField_impl.hpp
//---------------------------------------------------------------------------//
