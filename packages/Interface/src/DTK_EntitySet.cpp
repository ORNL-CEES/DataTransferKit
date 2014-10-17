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
 * \brief DTK_EntitySet.cpp
 * \author Stuart R. Slattery
 * \brief Geometric entity set interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_EntitySet.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_CommHelpers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
EntitySet::EntitySet()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
EntitySet::~EntitySet()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> > EntitySet::communicator() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return Teuchos::null;
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set. 
int EntitySet::physicalDimension() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}

//---------------------------------------------------------------------------//
// Get the local bounding box of entities of the set. Default implementation
// gathers the bounding boxes of local entities.
void EntitySet::localBoundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    EntityIterator entity_begin;
    EntityIterator entity_end;
    EntityIterator entity_it;
    EntityIterator dim_it;
    Teuchos::Tuple<double,6> entity_bounds;
    for ( int i = 0; i < 4; ++i )
    {
	dim_it = this->entityIterator(
	    static_cast<EntityType>(i),EntitySet::selectAll);
	entity_begin = dim_it.begin();
	entity_end = dim_it.end();
	for ( entity_it = entity_begin;
	      entity_it != entity_end;
	      ++entity_it )
	{
	    entity_it->boundingBox( entity_bounds );
	    for ( int n = 0; n < 3; ++n )
	    {
		bounds[n] = std::min( bounds[n], entity_bounds[n] );
		bounds[n+3] = std::max( bounds[n], entity_bounds[n] );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Get the global bounding box of entities of the set. Default implementation
// performs a parallel reduction using the local bounding boxes.
void EntitySet::globalBoundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    Teuchos::Tuple<double,6> local_bounds;
    this->localBoundingBox( local_bounds );
    Teuchos::reduceAll( *(this->communicator()), Teuchos::REDUCE_MIN, 3,
			&local_bounds[0], &bounds[0] ); 
    Teuchos::reduceAll( *(this->communicator()), Teuchos::REDUCE_MAX, 3,
			&local_bounds[3], &bounds[3] ); 
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void EntitySet::getEntity( const EntityId entity_id, Entity& entity ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate. 
EntityIterator EntitySet::entityIterator(
    const EntityType entity_type,
    const std::function<bool(Entity)>& predicate ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return EntityIterator();
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it. 
void EntitySet::getAdjacentEntities(
    const Entity& entity,
    const EntityType entity_type,
    Teuchos::Array<Entity>& adjacent_entities ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntitySet.cpp
//---------------------------------------------------------------------------//
