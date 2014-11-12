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
 * \brief DTK_STKMeshEntityImpl.cpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity implementation.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_STKMeshEntityImpl.hpp"
#include "DTK_DBC.hpp"

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
//brief Destructor.
STKMeshEntityImpl::~STKMeshEntityImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the entity type.
EntityType STKMeshEntityImpl::entityType() const
{
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    stk::mesh::EntityRank rank = 
	d_bulk_data->entity_rank(d_extra_data->d_stk_entity);
    
    EntityType entity_type = ENTITY_TYPE_NODE;
    switch( rank )
    {
	case stk::topology::NODE_RANK:
	    entity_type = ENTITY_TYPE_NODE;
	    break;
	case stk::topology::EDGE_RANK:
	    entity_type = ENTITY_TYPE_EDGE;
	    break;
	case stk::topology::FACE_RANK:
	    entity_type = ENTITY_TYPE_FACE;
	    break;
	case stk::topology::ELEM_RANK:
	    entity_type = ENTITY_TYPE_VOLUME;
	    break;
	default:
	    DTK_CHECK( ENTITY_TYPE_NODE == entity_type ||
		       ENTITY_TYPE_EDGE == entity_type ||
		       ENTITY_TYPE_FACE == entity_type ||
		       ENTITY_TYPE_VOLUME == entity_type );
	    break;
    }

    return entity_type;
}

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
    Teuchos::Array<stk::mesh::Entity> entity_nodes;
    stk::mesh::EntityRank rank = 
	d_bulk_data->entity_rank(d_extra_data->d_stk_entity);
    if ( stk::topology::NODE_RANK == rank )
    {
	entity_nodes.push_back( d_extra_data->d_stk_entity );
    }
    else
    {
	const stk::mesh::Entity* begin = 
	    d_bulk_data->begin_nodes( d_extra_data->d_stk_entity );
	const stk::mesh::Entity* end = 
	    d_bulk_data->end_nodes( d_extra_data->d_stk_entity );
	entity_nodes.assign( begin, end );
    }

    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    int space_dim = physicalDimension();
    switch( space_dim )
    {
	case 3:
	    getNodeBounds<stk::mesh::Cartesian3d>( entity_nodes, bounds );
	    break;
	case 2:
	    getNodeBounds<stk::mesh::Cartesian2d>( entity_nodes, bounds );
	    break;
	default:
	    DTK_CHECK( 2 == space_dim || 3 == space_dim );
	    break;
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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityImpl.cpp
//---------------------------------------------------------------------------//
