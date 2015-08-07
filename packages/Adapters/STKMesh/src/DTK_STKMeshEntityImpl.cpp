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
 * \brief STK mesh entity implementation.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_STKMeshEntityImpl.hpp"
#include "DTK_STKMeshHelpers.hpp"
#include "DTK_DBC.hpp"

#include <Intrepid_FieldContainer.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
STKMeshEntityImpl::STKMeshEntityImpl(
    const stk::mesh::Entity& stk_entity,
    const Teuchos::Ptr<stk::mesh::BulkData>& bulk_data )
    : d_extra_data( new STKMeshEntityExtraData(stk_entity) )
    , d_bulk_data( bulk_data )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId STKMeshEntityImpl::id() const
{ 
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    return Teuchos::as<EntityId>( 
	d_bulk_data->identifier(d_extra_data->d_stk_entity) );
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int STKMeshEntityImpl::ownerRank() const
{ 
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    return d_bulk_data->parallel_owner_rank( d_extra_data->d_stk_entity );
}

//---------------------------------------------------------------------------//
// Get the topological dimension of the entity.
int STKMeshEntityImpl::topologicalDimension() const
{
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    stk::mesh::EntityRank rank = 
	d_bulk_data->entity_rank(d_extra_data->d_stk_entity);
    return STKMeshHelpers::getTopologicalDimensionFromRank(
	rank, physicalDimension() );
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int STKMeshEntityImpl::physicalDimension() const
{ 
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    return d_bulk_data->mesh_meta_data().spatial_dimension();
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
void STKMeshEntityImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );

    Intrepid::FieldContainer<double> node_coords =
	STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(1,d_extra_data->d_stk_entity), 
	    *d_bulk_data );
    DTK_CHECK( node_coords.rank() == 3 );
    DTK_CHECK( node_coords.dimension(0) == 1 );

    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    int space_dim = node_coords.dimension(2);
    for ( int n = 0; n < node_coords.dimension(1); ++n )
    {
	for ( int d = 0; d < space_dim; ++d )
	{
	    bounds[d] = std::min( bounds[d], node_coords(0,n,d) );
	    bounds[d+3] = std::max( bounds[d+3], node_coords(0,n,d) );
	}
    }
    for ( int d = space_dim; d < 3; ++d )
    {
	bounds[d] = -max;
	bounds[d+3] = max;
    }

}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool STKMeshEntityImpl::inBlock( const int block_id ) const
{
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    const stk::mesh::PartVector& all_parts =
	d_bulk_data->mesh_meta_data().get_parts();
    stk::mesh::Bucket& entity_bucket =
	d_bulk_data->bucket( d_extra_data->d_stk_entity );
    for ( auto part_it = all_parts.begin(); 
	  part_it != all_parts.end();
	  ++part_it )
    {
	if ( Teuchos::as<int>((*part_it)->mesh_meta_data_ordinal()) == block_id )
	{
	    return entity_bucket.member( **part_it );
	}
    }
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool STKMeshEntityImpl::onBoundary( const int boundary_id ) const
{
    return inBlock( boundary_id );
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
Teuchos::RCP<EntityExtraData> STKMeshEntityImpl::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
void STKMeshEntityImpl::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    shards::CellTopology topo =
	STKMeshHelpers::getShardsTopology( d_extra_data->d_stk_entity,
					   *d_bulk_data );
    
    Intrepid::FieldContainer<double> node_coords =
	STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(1,d_extra_data->d_stk_entity), 
	    *d_bulk_data );
    int num_node = node_coords.dimension(1);
    int space_dim = node_coords.dimension(2);

    out << std::endl;
    out << "---" << std::endl;
    out << "STK Mesh Entity" << std::endl;
    out << "Id: " << id() << std::endl;
    out << "Owner rank: " << ownerRank() << std::endl;
    out << "Topology: " << topo;
    out << "Node coords: " << std::endl;
    for ( int n = 0; n < num_node; ++n )
    {
	out << "    node " << n << ": ";
	for ( int d = 0; d < space_dim; ++d )
	{
	    out << node_coords(0,n,d) << "  "; 
	}
	out << std::endl;
    }
    out << "---" << std::endl;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityImpl.cpp
//---------------------------------------------------------------------------//
