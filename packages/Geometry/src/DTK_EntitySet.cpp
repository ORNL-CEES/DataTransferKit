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
// Return a string indicating the derived entity set type.
std::string EntitySet::name() const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
    return "Not Implemented";
}

//---------------------------------------------------------------------------//
// Assign a parallel communicator to the entity set. This will only be done
// immediately after construct through the AbstractBuilder interface.
void EntitySet::assignCommunicator( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> > EntitySet::communicator() const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
    return Teuchos::null;
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate. 
AbstractIterator<GeometricEntity>
EntitySet::entityIterator(
    const EntityType entity_type,
    const std::function<bool(const GeometricEntity&)>& predicate ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
    return AbstractIterator<Entity>();
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void EntitySet::getEntity( const EntityId entity_id, 
			   Entity& entity ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Indicate that the entity set will be modified.
void EntitySet::startModification()
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Add an entity to the set.
void EntitySet::addEntity( const Entity& entity )
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Indicate that modification of the entity set is complete.
void EntitySet::endModification()
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set.
int EntitySet::physicalDimension() const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
    return -1;
}

//---------------------------------------------------------------------------//
// Get the local bounding box of entities of the set.
void EntitySet::localBoundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Get the global bounding box of entities of the set.
void EntitySet::globalBoundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntitySet.cpp
//---------------------------------------------------------------------------//
