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
// Create an empty clone of the EntitySet.
Teuchos::RCP<EntitySet> EntitySet::clone() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return Teuchos::null;
}

//---------------------------------------------------------------------------//
// Return a string indicating the derived entity set type.
std::string EntitySet::entitySetType() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return "Not Implemented";
}

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> > EntitySet::communicator() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return Teuchos::null;
}

//---------------------------------------------------------------------------//
// Get the local number of entities in the set of the given
std::size_t EntitySet::localNumberOfEntities( 
    const int parametric_dimension ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return 0;
}

//---------------------------------------------------------------------------//
// Get the global number of entities in the set of the given
std::size_t EntitySet::globalNumberOfEntities(
    const int parametric_dimension ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return 0;
}
    
//---------------------------------------------------------------------------//
// Get the identifiers for all local entities in the set of a given
void EntitySet::localEntityIds( 
    const int parametric_dimension,
    const Teuchos::ArrayView<EntityId>& entity_ids ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void EntitySet::getEntity( const EntityId entity_id, 
			   Teuchos::RCP<GeometricEntity>& entity ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Indicate that the entity set will be modified.
void EntitySet::startModification()
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Add an entity to the set.
void EntitySet::addEntity( const Teuchos::RCP<GeometricEntity>& entity )
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Indicate that modification of the entity set is complete.
void EntitySet::endModification()
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Return the largest spatial dimension of the entities in the set.
int EntitySet::spatialDimension() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}

//---------------------------------------------------------------------------//
// Get the local bounding box of entities of the set.
void EntitySet::localBoundingBox( Box& bounding_box ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Get the global bounding box of entities of the set.
void EntitySet::globalBoundingBox( Box& bounding_box ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntitySet.cpp
//---------------------------------------------------------------------------//
