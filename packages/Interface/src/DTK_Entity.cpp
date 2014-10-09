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
 * \brief DTK_Entity.cpp
 * \author Stuart R. Slattery
 * \brief Geometric entity interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_Entity.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
Entity::Entity()
{ /* ... */ }

//---------------------------------------------------------------------------//
//brief Destructor.
Entity::~Entity()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return a string indicating the derived entity type.
std::string Entity::name() const
{
    return b_entity_impl->name();
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId Entity::id() const
{ 
    return b_entity_impl->id();
}
    
//---------------------------------------------------------------------------//
// Get the entity type.
EntityType Entity::entityType() const
{
    return b_entity_impl->entityType();
}

//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int Entity::ownerRank() const
{ 
    return b_entity_impl->ownerRank();
}
//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int Entity::physicalDimension() const
{ 
    return b_entity_impl->physicalDimension();
}

//---------------------------------------------------------------------------//
// Return the axis-aligned bounding box around the entity.
void Entity::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{ 
    b_entity_impl->boundingBox( bounds );
}

//---------------------------------------------------------------------------//
// Return a string indicating the derived object type.
std::string Entity::objectType() const
{
    return b_entity_impl->name();
}

//---------------------------------------------------------------------------//
// Serialize the entity into a buffer.
void Entity::serialize( 
    const Teuchos::ArrayView<char>& buffer ) const
{
    b_entity_impl->serialize( buffer );
}

//---------------------------------------------------------------------------//
// Deserialize an entity from a buffer.
void Entity::deserialize( 
    const Teuchos::ArrayView<const char>& buffer )
{
    b_entity_impl->deserialize( buffer );
}

//---------------------------------------------------------------------------//
// Check whether the underlying implementation is available.
bool Entity::isEntityImplNonnull() const
{
    return Teuchos::nonnull( b_entity_impl );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Entity.cpp
//---------------------------------------------------------------------------//
