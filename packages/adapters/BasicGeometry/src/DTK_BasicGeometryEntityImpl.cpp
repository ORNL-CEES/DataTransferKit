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
 * \file DTK_BasicGeometryEntity.cpp
 * \author Stuart R. Slattery
 * \brief BasicGeometryEntity definition
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_DBC.hpp"
#include "DTK_BasicGeometryEntityImpl.hpp"
#include "DTK_BasicGeometryExtraData.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
BasicGeometryEntityImpl::BasicGeometryEntityImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
BasicGeometryEntityImpl::~BasicGeometryEntityImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the entity type.
EntityType BasicGeometryEntityImpl::entityType() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return static_cast<EntityType>(-1);
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId BasicGeometryEntityImpl::id() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int BasicGeometryEntityImpl::ownerRank() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}
//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int BasicGeometryEntityImpl::physicalDimension() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
void BasicGeometryEntityImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the surface of the set.
bool BasicGeometryEntityImpl::onSurface() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool BasicGeometryEntityImpl::inBlock( const int block_id ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool BasicGeometryEntityImpl::onBoundary( const int boundary_id ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
Teuchos::RCP<EntityExtraData> BasicGeometryEntityImpl::extraData() const
{
    return Teuchos::rcp( new BasicGeometryExtraData(this) );
}

//---------------------------------------------------------------------------//
// Compute the measure of the entity.
double BasicGeometryEntityImpl::measure() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1.0;
}

//---------------------------------------------------------------------------//
// Get the centroid of the entity.
void BasicGeometryEntityImpl::centroid( 
    const Teuchos::ArrayView<double>& centroid ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Safeguard the reverse map.
bool BasicGeometryEntityImpl::isSafeToMapToReferenceFrame(
    const Teuchos::ArrayView<const double>& point ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the
bool BasicGeometryEntityImpl::mapToReferenceFrame( 
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of
bool BasicGeometryEntityImpl::checkPointInclusion( 
    const double tolerance,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void BasicGeometryEntityImpl::mapToPhysicalFrame( 
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicGeometryEntityImpl.cpp
//---------------------------------------------------------------------------//

