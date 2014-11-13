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
 * \brief DTK_STKMeshHelpers.cpp
 * \author Stuart R. Slattery
 * \brief STK mesh helpers.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshHelpers.hpp"
#include "DTK_DBC.hpp"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Given a DTK EntityType, get the STK entity rank.
stk::mesh::EntityRank STKMeshHelpers::getRankFromType( 
    const EntityType dtk_type, const int space_dim )
{
    stk::mesh::EntityRank stk_rank = stk::topology::INVALID_RANK;

    switch( space_dim )
    {
	case 3:
	    switch( dtk_type )
	    {
		case ENTITY_TYPE_NODE:
		    stk_rank = stk::topology::NODE_RANK;
		    break;
		case ENTITY_TYPE_EDGE:
		    stk_rank = stk::topology::EDGE_RANK;
		    break;
		case ENTITY_TYPE_FACE:
		    stk_rank = stk::topology::FACE_RANK;
		    break;
		case ENTITY_TYPE_VOLUME:
		    stk_rank = stk::topology::ELEM_RANK;
		    break;
		default:
		    DTK_CHECK( ENTITY_TYPE_NODE == dtk_type ||
			       ENTITY_TYPE_EDGE == dtk_type ||
			       ENTITY_TYPE_FACE == dtk_type ||
			       ENTITY_TYPE_VOLUME == dtk_type );
		    break;
	    }
	    break;

	case 2:
	    switch( dtk_type )
	    {
		case ENTITY_TYPE_NODE:
		    stk_rank = stk::topology::NODE_RANK;
		    break;
		case ENTITY_TYPE_EDGE:
		    stk_rank = stk::topology::EDGE_RANK;
		    break;
		case ENTITY_TYPE_FACE:
		    stk_rank = stk::topology::ELEM_RANK;
		    break;
		default:
		    DTK_CHECK( ENTITY_TYPE_NODE == dtk_type ||
			       ENTITY_TYPE_EDGE == dtk_type ||
			       ENTITY_TYPE_FACE == dtk_type );
		    break;
	    }
	    break;

	default:
	    stk_rank = stk::topology::NODE_RANK;
    }

    return stk_rank;
}

//---------------------------------------------------------------------------//
// Given a STK entity stk_rank, get the DTK entity type.
EntityType STKMeshHelpers::getTypeFromRank(
    const stk::mesh::EntityRank stk_rank )
{
    EntityType dtk_type = ENTITY_TYPE_NODE;
    switch( stk_rank )
    {
	case stk::topology::NODE_RANK:
	    dtk_type = ENTITY_TYPE_NODE;
	    break;
	case stk::topology::EDGE_RANK:
	    dtk_type = ENTITY_TYPE_EDGE;
	    break;
	case stk::topology::FACE_RANK:
	    dtk_type = ENTITY_TYPE_FACE;
	    break;
	case stk::topology::ELEM_RANK:
	    dtk_type = ENTITY_TYPE_VOLUME;
	    break;
	default:
	    DTK_CHECK( ENTITY_TYPE_NODE == dtk_type ||
		       ENTITY_TYPE_EDGE == dtk_type ||
		       ENTITY_TYPE_FACE == dtk_type ||
		       ENTITY_TYPE_VOLUME == dtk_type );
	    break;
    }

    return dtk_type;
}

//---------------------------------------------------------------------------//
// Given a DTK entity, return the corresponding STK entity key.
stk::mesh::EntityKey 
STKMeshHelpers::getKeyFromEntity( const Entity dtk_entity )
{
    return stk::mesh::EntityKey( 
	getRankFromType(dtk_entity.entityType(),dtk_entity.physicalDimension()),
	dtk_entity.id() );
}

//---------------------------------------------------------------------------//
// Given a STK entity, return the coordinates of its nodes in a field
// container ordered by canonical node order (N,D).
Intrepid::FieldContainer<double> 
STKMeshHelpers::getEntityNodeCoordinates( 
    const Teuchos::Array<stk::mesh::Entity>& stk_entities,
    const stk::mesh::BulkData& bulk_data )
{
    int space_dim = bulk_data.mesh_meta_data().spatial_dimension();
    switch( space_dim )
    {
	case 3:
	    return
		STKMeshHelpers::extractEntityNodeCoordinates<
		    stk::mesh::Cartesian3d>( stk_entities, bulk_data, space_dim );
	    break;
	case 2:
	    return
		STKMeshHelpers::extractEntityNodeCoordinates<
		    stk::mesh::Cartesian2d>( stk_entities, bulk_data, space_dim );
	    break;
	default:
	    DTK_CHECK( 2 == space_dim || 3 == space_dim );
	    break;
    }
    return Intrepid::FieldContainer<double>();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshHelpers.cpp
//---------------------------------------------------------------------------//
