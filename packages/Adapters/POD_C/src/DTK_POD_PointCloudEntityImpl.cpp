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
 * \brief DTK_STKMeshEntityImpl.cpp
 * \author Stuart R. Slattery
 * \brief Point cloud entity implementation.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_POD_PointCloudEntityImpl.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
POD_PointCloudEntityImpl::POD_PointCloudEntityImpl(
    const double* cloud_coords,
    const unsigned num_points,
    const int space_dim,
    const DTK_Data_layout layout,
    const EntityId global_id,                          
    const int local_id,
    const int owner_rank )
    : d_cloud_coords( cloud_coords )
    , d_offsets( space_dim, -1 )
    , d_global_id( global_id )
    , d_owner_rank( owner_rank )
{
    DTK_REQUIRE( DTK_INTERLEAVED == layout ||
                 DTK_BLOCKED == layout );
    
    // Calculate the offsets into the coordinates array.
    for ( int d = 0; d < space_dim; ++d )
    {
        d_offsets[d] = ( DTK_INTERLEAVED == layout )
                       ? space_dim*local_id + d
                       : d*num_points + local_id;
    }
}

//---------------------------------------------------------------------------//
// Get the coordinates of the point in a given dimension.
double POD_PointCloudEntityImpl::coords( const int dim ) const
{
    DTK_REQUIRE( dim < physicalDimension() );
    return d_cloud_coords[ d_offsets[dim] ];
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId POD_PointCloudEntityImpl::id() const
{ 
    return d_global_id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int POD_PointCloudEntityImpl::ownerRank() const
{ 
    return d_owner_rank;
}

//---------------------------------------------------------------------------//
// Get the topological dimension of the entity.
int POD_PointCloudEntityImpl::topologicalDimension() const
{
    return 0;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int POD_PointCloudEntityImpl::physicalDimension() const
{ 
    return d_offsets.size();
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
void POD_PointCloudEntityImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    for ( int d = 0; d < physicalDimension(); ++d )
    {
        bounds[d] = d_cloud_coords[ d_offsets[d] ];
        bounds[d+3] = bounds[d];
    }

    for ( int d = physicalDimension(); d < 3; ++d )
    {
        bounds[d] = 0.0;
        bounds[d+3] = bounds[d];
    }
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool POD_PointCloudEntityImpl::inBlock( const int block_id ) const
{
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool POD_PointCloudEntityImpl::onBoundary( const int boundary_id ) const
{
    return true;
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
Teuchos::RCP<EntityExtraData> POD_PointCloudEntityImpl::extraData() const
{
    return Teuchos::null;
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
void POD_PointCloudEntityImpl::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{

    out << std::endl;
    out << "---" << std::endl;
    out << "POD Point Cloud Entity" << std::endl;

    out << "Id: " << id() << std::endl;

    out << "Owner rank: " << ownerRank() << std::endl;

    out << "Point coords: ";
    for ( int d = 0; d < physicalDimension(); ++d )
    {
        out << d_cloud_coords[d_offsets[d]] << " ";
    }
    out << std::endl;        
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_POD_PointCloudEntityImpl.cpp
//---------------------------------------------------------------------------//
