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
 * \brief DTK_STKMeshEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief Geometric entity set interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshEntitySet.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
STKMeshEntitySet::STKMeshEntitySet( 
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data )
    : d_bulk_data( bulk_data )
{ 
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
}

//---------------------------------------------------------------------------//
// Destructor.
STKMeshEntitySet::~STKMeshEntitySet()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> > STKMeshEntitySet::communicator() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return Teuchos::null;
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set. 
int STKMeshEntitySet::physicalDimension() const
{
    return d_bulk_data->mesh_meta_data().spatial_dimension();
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void STKMeshEntitySet::getEntity( const EntityType entity_type,
				  const EntityId entity_id, 
				  Entity& entity ) const
{

}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate. 
EntityIterator STKMeshEntitySet::entityIterator(
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
void STKMeshEntitySet::getAdjacentEntities(
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
// end DTK_STKMeshEntitySet.cpp
//---------------------------------------------------------------------------//
