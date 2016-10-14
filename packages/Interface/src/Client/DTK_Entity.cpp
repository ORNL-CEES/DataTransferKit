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
 * \brief DTK_Entity.cpp
 * \author Stuart R. Slattery
 * \brief Geometric entity interface.
 */
//---------------------------------------------------------------------------//

#include <sstream>

#include "DTK_DBC.hpp"
#include "DTK_Entity.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
Entity::Entity() { /* ... */}

//---------------------------------------------------------------------------//
// Copy constructor.
Entity::Entity( const Entity &rhs ) { b_entity_impl = rhs.b_entity_impl; }

//---------------------------------------------------------------------------//
// Copy assignment operator.
Entity &Entity::operator=( const Entity &rhs )
{
    b_entity_impl = rhs.b_entity_impl;
    return *this;
}

//---------------------------------------------------------------------------//
// Move constructor.
Entity::Entity( Entity &&rhs ) { b_entity_impl = rhs.b_entity_impl; }

//---------------------------------------------------------------------------//
// Move assignment operator.
Entity &Entity::operator=( Entity &&rhs )
{
    b_entity_impl = rhs.b_entity_impl;
    return *this;
}

//---------------------------------------------------------------------------//
// brief Destructor.
Entity::~Entity() { /* ... */}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId Entity::id() const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->id();
}

//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int Entity::ownerRank() const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->ownerRank();
}

//---------------------------------------------------------------------------//
// Return the topological dimension of the entity.
int Entity::topologicalDimension() const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->topologicalDimension();
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int Entity::physicalDimension() const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->physicalDimension();
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
void Entity::boundingBox( Teuchos::Tuple<double, 6> &bounds ) const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    b_entity_impl->boundingBox( bounds );
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool Entity::inBlock( const int block_id ) const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->inBlock( block_id );
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool Entity::onBoundary( const int boundary_id ) const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->onBoundary( boundary_id );
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
Teuchos::RCP<EntityExtraData> Entity::extraData() const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->extraData();
}

//---------------------------------------------------------------------------//
// Provide a one line description of the object.
std::string Entity::description() const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    std::stringstream d;
    d << " Id = " << id() << ", OwnerRank = " << ownerRank()
      << ", TopologicalDimension = " << topologicalDimension()
      << ", PhysicalDimension = " << physicalDimension();
    return b_entity_impl->description() + d.str();
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
void Entity::describe( Teuchos::FancyOStream &out,
                       const Teuchos::EVerbosityLevel verb_level ) const
{
    DTK_REQUIRE( Teuchos::nonnull( b_entity_impl ) );
    return b_entity_impl->describe( out, verb_level );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Entity.cpp
//---------------------------------------------------------------------------//
