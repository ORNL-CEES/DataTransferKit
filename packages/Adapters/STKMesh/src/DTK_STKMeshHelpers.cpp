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
 * \brief DTK_STKMeshHelpers.cpp
 * \author Stuart R. Slattery
 * \brief STK mesh helpers.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshHelpers.hpp"
#include "DTK_STKMeshEntityExtraData.hpp"
#include "DTK_DBC.hpp"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_topology/topology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Given a DTK entity, extract the STK entity.
const stk::mesh::Entity& STKMeshHelpers::extractEntity( const Entity dtk_entity )
{
    return Teuchos::rcp_dynamic_cast<STKMeshEntityExtraData>(
        dtk_entity.extraData() )->d_stk_entity;
}

//---------------------------------------------------------------------------//
// Given a topological dimension, get the STK entity rank.
stk::mesh::EntityRank STKMeshHelpers::getRankFromTopologicalDimension(
    const int topo_dim, const int space_dim )
{
    stk::mesh::EntityRank stk_rank = stk::topology::INVALID_RANK;

    switch( space_dim )
    {
        case 3:
            switch( topo_dim )
            {
                case 0:
                    stk_rank = stk::topology::NODE_RANK;
                    break;
                case 1:
                    stk_rank = stk::topology::EDGE_RANK;
                    break;
                case 2:
                    stk_rank = stk::topology::FACE_RANK;
                    break;
                case 3:
                    stk_rank = stk::topology::ELEM_RANK;
                    break;
                default:
                    DTK_CHECK( 0 == topo_dim ||
                               1 == topo_dim ||
                               2 == topo_dim ||
                               3 == topo_dim );
                    break;
            }
            break;

        case 2:
            switch( topo_dim )
            {
                case 0:
                    stk_rank = stk::topology::NODE_RANK;
                    break;
                case 1:
                    stk_rank = stk::topology::EDGE_RANK;
                    break;
                case 2:
                    stk_rank = stk::topology::ELEM_RANK;
                    break;
                default:
                    DTK_CHECK( 0 == topo_dim ||
                               1 == topo_dim ||
                               2 == topo_dim );
                    break;
            }
            break;

        default:
            stk_rank = stk::topology::NODE_RANK;
    }

    return stk_rank;
}

//---------------------------------------------------------------------------//
// Given a STK entity stk_rank, get the topological dimension.
int STKMeshHelpers::getTopologicalDimensionFromRank(
    const stk::mesh::EntityRank stk_rank, const int space_dim )
{
    int topo_dim = 0;

    switch( space_dim )
    {
        case 3:
            switch( stk_rank )
            {
                case stk::topology::NODE_RANK:
                    topo_dim = 0;
                    break;
                case stk::topology::EDGE_RANK:
                    topo_dim = 1;
                    break;
                case stk::topology::FACE_RANK:
                    topo_dim = 2;
                    break;
                case stk::topology::ELEM_RANK:
                    topo_dim = 3;
                    break;
                default:
                    DTK_CHECK( stk::topology::NODE_RANK == stk_rank ||
                               stk::topology::EDGE_RANK == stk_rank ||
                               stk::topology::FACE_RANK == stk_rank ||
                               stk::topology::ELEM_RANK == stk_rank );
                    break;
            }
            break;

        case 2:
            switch( stk_rank )
            {
                case stk::topology::NODE_RANK:
                    topo_dim = 0;
                    break;
                case stk::topology::EDGE_RANK:
                    topo_dim = 1;
                    break;
                case stk::topology::ELEM_RANK:
                    topo_dim = 2;
                    break;
                default:
                    DTK_CHECK( stk::topology::NODE_RANK == stk_rank ||
                               stk::topology::EDGE_RANK == stk_rank ||
                               stk::topology::ELEM_RANK == stk_rank );
                    break;
            }
            break;

        default:
            DTK_CHECK( 3 == space_dim || 2 == space_dim );
            break;
    }

    return topo_dim;
}

//---------------------------------------------------------------------------//
// Given a DTK entity, return the corresponding STK entity key.
stk::mesh::EntityKey
STKMeshHelpers::getKeyFromEntity( const Entity dtk_entity )
{
    return stk::mesh::EntityKey(
        getRankFromTopologicalDimension(
            dtk_entity.topologicalDimension(),dtk_entity.physicalDimension()),
        dtk_entity.id() );
}

//---------------------------------------------------------------------------//
// Given a STK entity, return its shards topology.
shards::CellTopology
STKMeshHelpers::getShardsTopology( const stk::mesh::Entity stk_entity,
                                   const stk::mesh::BulkData& bulk_data )
{
    return stk::mesh::get_cell_topology(
        bulk_data.bucket(stk_entity).topology() );
}

//---------------------------------------------------------------------------//
// Given a set of STK entities, return the coordinates of its nodes in a field
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
