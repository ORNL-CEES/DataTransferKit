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
#include <algorithm>

#include "DTK_DBC.hpp"
#include "DTK_PointImpl.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
PointImpl::PointImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Array constructor.
PointImpl::PointImpl( const EntityId global_id, 
		      const int owner_rank,
		      const Teuchos::Array<double>& coordinates,
		      const bool on_surface,
		      const Teuchos::Array<int>& block_ids,
		      const Teuchos::Array<int>& boundary_ids )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_on_surface( on_surface )
    , d_block_ids( block_ids )
    , d_boundary_ids( boundary_ids )
    , d_coordinates( coordinates )
{
    std::sort( d_block_ids.begin(), d_block_ids.end() );
    std::sort( d_boundary_ids.begin(), d_boundary_ids.end() );
}

//---------------------------------------------------------------------------//
// Destructor.
PointImpl::~PointImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the coordinates of the point.
void PointImpl::getCoordinates( 
    const Teuchos::ArrayView<double>& coordinates ) const
{ 
    coordinates.assign( d_coordinates );
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
// Return the axis-aligned bounding box around the entity.
void PointImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    Teuchos::Array<double> coordinates( this->physicalDimension() );
    this->getCoordinates( coordinates );
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
// Determine if an entity is on the surface of the set.
bool PointImpl::onSurface() const
{
    return d_on_surface;
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool PointImpl::inBlock( const int block_id ) const
{
    return std::binary_search( 
	d_block_ids.begin(), d_block_ids.end(), block_id );
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool PointImpl::onBoundary( const int boundary_id ) const
{
    return std::binary_search( 
	d_boundary_ids.begin(), d_boundary_ids.end(), boundary_id );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the measure of the point.
 *
 * \return Return the measure of the point.
 */
double PointImpl::measure() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the point.
 *
 * \return The centroid coordinates.
 */
void PointImpl::centroid( const Teuchos::ArrayView<double>& centroid ) const
{
    this->getCoordinates( centroid );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Safeguard the reverse map.
 */
bool PointImpl::isSafeToMapToReferenceFrame(
    const Teuchos::ArrayView<const double>& point ) const
{
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a point to the reference space of an entity. Return the
 */
bool PointImpl::mapToReferenceFrame( 
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    reference_point.assign( point );
    return false;
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Determine if a reference point is in the parameterized space of
 * an entity. This is true only if the points are separated by a distance less
 * than the tolerance.
 */
bool PointImpl::checkPointInclusion( 
    const double tolerance,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    int space_dim = this->physicalDimension();
    Teuchos::Array<double> coords( space_dim );
    this->getCoordinates( coords() );
    double distance = 0.0;
    double local_dist = 0.0;
    for ( int d = 0; d < space_dim; ++d )
    {
	local_dist = coords[d] - reference_point[d];
	distance += local_dist*local_dist;
    }
    return ( distance < tolerance*tolerance );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a reference point to the physical space of an entity.
 */
void PointImpl::mapToPhysicalFrame( 
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    point.assign( reference_point );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_PointImpl.cpp
//---------------------------------------------------------------------------//

