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

#include "DTK_ReferenceHexImpl.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
namespace UnitTest
{
//---------------------------------------------------------------------------//
// Constructor.
ReferenceHexImpl::ReferenceHexImpl(
    const int id,
    const int owner_rank,
    const Teuchos::Array<DataTransferKit::Entity>& nodes )
    : d_extra_data( new ReferenceHexExtraData() )
{
    DTK_REQUIRE( 8 == nodes.size() );
    
    // Get the element id and owner rank.
    d_extra_data->id = id;
    d_extra_data->owner_rank = owner_rank;

    // Get the element node ids.
    d_extra_data->node_ids.resize( 8 );
    for ( int n = 0; n < 8; ++n )
    {
        d_extra_data->node_ids[n] = nodes[n].id();
    }

    // Get the element node coordinates.
    Teuchos::Tuple<double,6> node_bounds;
    auto& coords = d_extra_data->node_coords;
    coords.resize( 1, 8, 3 );
    for ( int n = 0; n < 8; ++n )
    {
        nodes[n].boundingBox( node_bounds );
        for ( int d = 0; d < 3; ++d )
        {
            coords(0,n,d) = node_bounds[d];
        }
    }
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
DataTransferKit::EntityId ReferenceHexImpl::id() const
{ 
    return d_extra_data->id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int ReferenceHexImpl::ownerRank() const
{ 
    return d_extra_data->owner_rank;
}

//---------------------------------------------------------------------------//
// Get the topological dimension of the entity.
int ReferenceHexImpl::topologicalDimension() const
{
    return 3;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int ReferenceHexImpl::physicalDimension() const
{ 
    return 3;
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
void ReferenceHexImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    for ( int n = 0; n < 8; ++n )
    {
	for ( int d = 0; d < 3; ++d )
	{
	    bounds[d] = std::min( bounds[d],
                                  d_extra_data->node_coords(0,n,d) );
	    bounds[d+3] = std::max( bounds[d+3],
                                    d_extra_data->node_coords(0,n,d) );
	}
    }
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool ReferenceHexImpl::inBlock( const int block_id ) const
{
    // no blocks for now.
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool ReferenceHexImpl::onBoundary( const int boundary_id ) const
{
    // no boundaries for now.
    return false;
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
Teuchos::RCP<DataTransferKit::EntityExtraData>
ReferenceHexImpl::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
void ReferenceHexImpl::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    out << std::endl;
    out << "---" << std::endl;
    out << "Reference Hex Entity" << std::endl;
    out << "Id: " << id() << std::endl;
    out << "Owner rank: " << ownerRank() << std::endl;
    out << "Node ids and coords: " << std::endl;
    for ( int n = 0; n < 8; ++n )
    {
	out << "    node " << n << ", id "
            << d_extra_data->node_ids[n] << ": ";
	for ( int d = 0; d < 3; ++d )
	{
	    out << d_extra_data->node_coords(0,n,d) << "  "; 
	}
	out << std::endl;
    }
    out << "---" << std::endl;
}

//---------------------------------------------------------------------------//

} // end namespace UnitTest
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ReferenceHexImpl.cpp
//---------------------------------------------------------------------------//
