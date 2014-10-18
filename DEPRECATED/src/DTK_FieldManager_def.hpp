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
 * \file DTK_FieldManager_def.hpp
 * \author Stuart R. Slattery
 * \brief Field manager definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDMANAGER_DEF_HPP
#define DTK_FIELDMANAGER_DEF_HPP

#include "DTK_DBC.hpp"
#include "DTK_FieldTools.hpp"
#include "DataTransferKit_config.hpp"

#include <Teuchos_CommHelpers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor. If Design-By-Contract is enabled, the constructor will
 * validate the field description to the domain model. This requires a few
 * global communications.
 *
 * \param field The field that this object is managing. This field must have
 * FieldTraits.
 * 
 * \param comm The communicator over which the field is defined.
 */
template<class Field>
FieldManager<Field>::FieldManager( const RCP_Field& field, 
				   const RCP_Comm& comm )
    : d_field( field )
    , d_comm( comm )
{
    // If we're checking with Design-by-Contract, validate the field to the
    // domain model.
#if HAVE_DTK_DBC
    validate();
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Field>
FieldManager<Field>::~FieldManager()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Validate the field to the domain model.
 */
template<class Field>
void FieldManager<Field>::validate()
{
    // Check that the field dimension is the same on every node.
    Teuchos::Array<int> local_dims( d_comm->getSize(), 0 );
    Teuchos::Array<int> local_dims_copy( d_comm->getSize(), 0 );
    local_dims[ d_comm->getRank() ] = FT::dim( *d_field );
    Teuchos::reduceAll<int,int>( *d_comm, Teuchos::REDUCE_SUM,
				 local_dims.size(),
				 &local_dims[0], &local_dims_copy[0] ); 
    Teuchos::Array<int>::iterator unique_bound;
    std::sort( local_dims_copy.begin(), local_dims_copy.end() );
    unique_bound = std::unique( local_dims_copy.begin(), local_dims_copy.end() );
    int unique_dim = std::distance( local_dims_copy.begin(), unique_bound );
    DTK_REQUIRE( 1 == unique_dim );
    local_dims_copy.clear();

    // Check that the data dimension is the same as the field dimension.
    typename FT::size_type num_data = std::distance( FT::begin( *d_field ), 
						     FT::end( *d_field ) );
    DTK_REQUIRE( num_data == FT::size( *d_field ) );
    if ( !FT::empty( *d_field ) )
    {
	DTK_REQUIRE( num_data / FieldTools<Field>::dimSize( *d_field ) 
			  == Teuchos::as<typename FT::size_type>(
			      FT::dim(*d_field)) );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_FIELDMANAGER_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldManager_def.hpp
//---------------------------------------------------------------------------//

