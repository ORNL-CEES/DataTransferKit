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
void EntityPredicates::setBlockId( const int block_id ) 
{ d_block_id = block_id; }
bool EntityPredicates::inBlock( Entity entity ) 
{ return entity.inBlock(d_block_id); }

//---------------------------------------------------------------------------//
// Boundary predicate.
void EntityPredicates::setBoundaryId( const int boundary_id ) 
{ d_boundary_id = boundary_id; }
bool EntityPredicates::onBoundary( Entity entity ) 
{ return entity.onBoundary(d_boundary_id); }

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
{ return (entity.entityType == d_entity_type); }

//---------------------------------------------------------------------------//
// Current block to check.
int EntityPredicates::d_block_id = -1;

//---------------------------------------------------------------------------//
// Current boundary to check.
int EntityPredicates::d_boundary_id = -1;

//---------------------------------------------------------------------------//
// Current owner rank to check.
int EntityPredicates::d_owner_rank = -1;

//---------------------------------------------------------------------------//
// Current entity type to check.
EntityType EntityPredicates::d_entity_type = ENTITY_TYPE_NODE

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityPredicates.cpp
//---------------------------------------------------------------------------//
