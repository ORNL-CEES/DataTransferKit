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
 * \brief DTK_EntityPredicates.cpp
 * \author Stuart R. Slattery
 * \brief Basic entity predicates.
 */
//---------------------------------------------------------------------------//

#include "DTK_EntityPredicates.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
//! Surface predicate.
bool EntityPredicates::onSurface( Entity entity ) 
{ return entity.onSurface(); }

//---------------------------------------------------------------------------//
// Block predicate.
void EntityPredicates::setBlockIds( const Teuchos::Array<int>& block_ids ) 
{ d_block_ids = block_ids; }
bool EntityPredicates::inBlocks( Entity entity ) 
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
// Boundary predicate.
void EntityPredicates::setBoundaryIds( const Teuchos::Array<int>& boundary_ids ) 
{ d_boundary_ids = boundary_ids; }
bool EntityPredicates::onBoundaries( Entity entity ) 
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
// Owner rank predicate.
void EntityPredicates::setOwnerRank( const int owner_rank )
{ d_owner_rank = owner_rank; }
bool EntityPredicates::hasOwner( Entity entity )
{ return (entity.ownerRank() == d_owner_rank); }

//---------------------------------------------------------------------------//
// Entity type predicate.
void EntityPredicates::setEntityType( const EntityType entity_type )
{ d_entity_type = entity_type; }
bool EntityPredicates::isEntityType( Entity entity )
{ return (entity.entityType() == d_entity_type); }

//---------------------------------------------------------------------------//
// Current blocks to check.
Teuchos::Array<int> EntityPredicates::d_block_ids = Teuchos::Array<int>();

//---------------------------------------------------------------------------//
// Current boundary to check.
Teuchos::Array<int> EntityPredicates::d_boundary_ids = Teuchos::Array<int>();

//---------------------------------------------------------------------------//
// Current owner rank to check.
int EntityPredicates::d_owner_rank = -1;

//---------------------------------------------------------------------------//
// Current entity type to check.
EntityType EntityPredicates::d_entity_type = ENTITY_TYPE_NODE;

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityPredicates.cpp
//---------------------------------------------------------------------------//
