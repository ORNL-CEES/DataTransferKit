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
 * \file DTK_Point.cpp
 * \author Stuart R. Slattery
 * \brief Point definition.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
PointImpl::PointImpl()
    : d_global_id( dtk_invalid_entity_id )
    , d_owner_rank( -1 )
    , d_coordinates( 0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Array constructor.
PointImpl::PointImpl( const EntityId global_id, 
			   const int owner_rank,
			   const Teuchos::Array<double>& coordinates )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_coordinates( coordinates )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
PointImpl::~PointImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the coordinates of the point.
void PointImpl::getCoordinates( 
    Teuchos::ArrayView<const double>& coordinates ) const
{ 
    coordinates = d_coordinates(); 
}

//---------------------------------------------------------------------------//
// Get the entity type.
EntityType PointImpl::entityType() const
{
    return ENTITY_TYPE_NODE;
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId PointImpl::id() const
{
    return d_global_id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int PointImpl::ownerRank() const
{
    return d_owner_rank;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int PointImpl::physicalDimension() const
{
    return d_coordinates.size();
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric
double PointImpl::measure() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void PointImpl::centroid( Teuchos::ArrayView<const double>& centroid ) const
{
    getCoordinates( centroid );
}

//---------------------------------------------------------------------------//
// Return the axis-aligned bounding box around the entity. 1D specialization.
void PointImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    Teuchos::ArrayView<const double> coordinates;
    getCoordinates( coordinates );
    double max = std::numeric_limits<double>::max();
    int space_dim = coordinates.size();
    if ( 1 == space_dim )
    {
	bounds = Teuchos::tuple( coordinates[0], -max, -max,
				 coordinates[0], max, max );
    }
    else if ( 2 == space_dim )
    {
	bounds = Teuchos::tuple( coordinates[0], coordinates[1], -max,
				 coordinates[0], coordinates[1], max );
    }
    else if ( 3 == space_dim )
    {
	bounds = Teuchos::tuple( coordinates[0], coordinates[1], coordinates[2], 
				 coordinates[0], coordinates[1], coordinates[2] );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_PointImpl.cpp
//---------------------------------------------------------------------------//

