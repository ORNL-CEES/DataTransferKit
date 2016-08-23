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
 * \brief DTK_ReferenceHexImpl.cpp
 * \author Stuart R. Slattery
 * \brief Reference hex mesh entity implementation.
 */
//---------------------------------------------------------------------------//

#include "DTK_ReferenceNodeImpl.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
namespace UnitTest
{
//---------------------------------------------------------------------------//
// Constructor.
ReferenceNodeImpl::ReferenceNodeImpl(
    const int id,
    const int owner_rank,
    const double x,
    const double y,
    const double z )    
    : d_extra_data( new ReferenceNodeExtraData() )
{
    d_extra_data->id = id;
    d_extra_data->owner_rank = owner_rank;
    d_extra_data->node_coords.resize( 3 );
    d_extra_data->node_coords[0] = x;
    d_extra_data->node_coords[1] = y;    
    d_extra_data->node_coords[2] = z;    
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
DataTransferKit::EntityId ReferenceNodeImpl::id() const
{ 
    return d_extra_data->id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int ReferenceNodeImpl::ownerRank() const
{ 
    return d_extra_data->owner_rank;
}

//---------------------------------------------------------------------------//
// Get the topological dimension of the entity.
int ReferenceNodeImpl::topologicalDimension() const
{
    return 0;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int ReferenceNodeImpl::physicalDimension() const
{ 
    return 3;
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
void ReferenceNodeImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( d_extra_data->node_coords[0],
                             d_extra_data->node_coords[1],
                             d_extra_data->node_coords[2],
                             d_extra_data->node_coords[0],
                             d_extra_data->node_coords[1],
                             d_extra_data->node_coords[2] );
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool ReferenceNodeImpl::inBlock( const int block_id ) const
{
    // no blocks for now.
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool ReferenceNodeImpl::onBoundary( const int boundary_id ) const
{
    // no boundaries for now.
    return false;
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
Teuchos::RCP<DataTransferKit::EntityExtraData>
ReferenceNodeImpl::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
void ReferenceNodeImpl::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    out << std::endl;
    out << "---" << std::endl;
    out << "Reference Node Entity" << std::endl;
    out << "Id: " << id() << std::endl;
    out << "Owner rank: " << ownerRank() << std::endl;
    out << "Node coords: " << std::endl;
    for ( int d = 0; d < 3; ++d )
    {
        out << d_extra_data->node_coords[d] << "  "; 
    }
    out << std::endl;
    out << "---" << std::endl;
}

//---------------------------------------------------------------------------//

} // end namespace UnitTest
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ReferenceNodeImpl.cpp
//---------------------------------------------------------------------------//
