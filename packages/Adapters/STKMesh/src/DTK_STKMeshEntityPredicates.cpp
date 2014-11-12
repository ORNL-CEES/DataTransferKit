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
 * \brief DTK_STKMeshEntityPredicates.cpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity predicates.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshEntityPredicates.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Part block predicate.
//---------------------------------------------------------------------------//
// Constructor.
PartBlockPredicate::PartBlockPredicate( 
    const Teuchos::Array<std::string>& part_names,
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data )
{
    for ( auto name_it = part_names.begin(); 
	  name_it != part_names.end();
	  ++name_it )
    {
	Part& part = bulk_data->mesh_meta_data.get_part( *name_it );
	d_block_ids.push_back( part.mesh_meta_data_ordinal() );
    }
}

//---------------------------------------------------------------------------//
// Functor.
bool PartBlockPredicate::operator()( Entity entity ) 
{ 
    Teuchos::Array<int>::const_iterator block_it;
    for ( block_it = d_block_ids.begin();
	  block_it != d_block_ids.end();
	  ++block_it )
    {
	if ( !entity.inBlock(*block_it) )
	{
	    return false;
	}
    }
    return true;
}

//---------------------------------------------------------------------------//
// Part boundary predicate.
//---------------------------------------------------------------------------//
// Constructor.
PartBoundaryPredicate::PartBoundaryPredicate( 
    const Teuchos::Array<std::string>& part_names,
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data )
{
    for ( auto name_it = part_names.begin(); 
	  name_it != part_names.end();
	  ++name_it )
    {
	Part& part = bulk_data->mesh_meta_data.get_part( *name_it );
	d_boundary_ids.push_back( part.mesh_meta_data_ordinal() );
    }
}

//---------------------------------------------------------------------------//
// Functor.
bool PartBoundaryPredicate::operator()( Entity entity ) 
{ 
    Teuchos::Array<int>::const_iterator boundary_it;
    for ( boundary_it = d_boundary_ids.begin();
	  boundary_it != d_boundary_ids.end();
	  ++boundary_it )
    {
	if ( !entity.onBoundary(*boundary_it) )
	{
	    return false;
	}
    }
    return true;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityPredicates.cpp
//---------------------------------------------------------------------------//
