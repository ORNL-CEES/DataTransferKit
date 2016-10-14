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
 * \brief DTK_MoabHelpers.cpp
 * \author Stuart R. Slattery
 * \brief Moab helper functions.
 */
//---------------------------------------------------------------------------//

#include "DTK_MoabHelpers.hpp"
#include "DTK_DBC.hpp"
#include "DTK_MoabEntityExtraData.hpp"

#include <MBTagConventions.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Given a DTK entity, extract the Moab entity.
moab::EntityHandle MoabHelpers::extractEntity( const Entity dtk_entity )
{
    return Teuchos::rcp_dynamic_cast<MoabEntityExtraData>(
               dtk_entity.extraData() )
        ->d_moab_entity;
}

//---------------------------------------------------------------------------//
// Get the global id of a list of entities.
void MoabHelpers::getGlobalIds( const moab::ParallelComm &mesh,
                                const moab::EntityHandle *entities,
                                const int num_entities, EntityId *entity_ids )
{
    // Get the global id tag.
    moab::Tag id_tag;
    DTK_CHECK_ERROR_CODE( mesh.get_moab()->tag_get_handle(
        GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag ) );

    // Extract the global ids. Note that if the Moab native global id
    // functionality is used as we are here then only ids of type 'int' are
    // available. That means no more than 3 billion entities may be uniquely
    // indexed in a mesh. For large problems we will either have to get Moab
    // updated to tag the mesh using a 'bit' type or construct global ids
    // ourselves.
    Teuchos::Array<int> id_tag_data( num_entities );
    DTK_CHECK_ERROR_CODE( mesh.get_moab()->tag_get_data(
        id_tag, entities, num_entities,
        static_cast<void *>( id_tag_data.getRawPtr() ) ) );

    // Copy the integer data into the array of type EntityId.
    std::copy( id_tag_data.begin(), id_tag_data.end(), entity_ids );
}

//---------------------------------------------------------------------------//
// Given a Moab EntityType, get the topological dimension.
int MoabHelpers::getTopologicalDimensionFromMoabType(
    const moab::EntityType moab_type )
{
    int topo_dim = 0;
    switch ( moab_type )
    {
    case moab::MBVERTEX:
        topo_dim = 0;
        break;
    case moab::MBEDGE:
        topo_dim = 1;
        break;
    case moab::MBTRI:
        topo_dim = 2;
        break;
    case moab::MBQUAD:
        topo_dim = 2;
        break;
    case moab::MBPOLYGON:
        topo_dim = 2;
        break;
    case moab::MBTET:
        topo_dim = 3;
        break;
    case moab::MBPYRAMID:
        topo_dim = 3;
        break;
    case moab::MBPRISM:
        topo_dim = 3;
        break;
    case moab::MBKNIFE:
        topo_dim = 3;
        break;
    case moab::MBHEX:
        topo_dim = 3;
        break;
    case moab::MBPOLYHEDRON:
        topo_dim = 3;
        break;
    case moab::MBENTITYSET:
        topo_dim = -1;
        break;
    default:
        topo_dim = -1;
        break;
    }
    return topo_dim;
}

//---------------------------------------------------------------------------//
// Given a Moab EntityType, get the topological dimension.
std::string MoabHelpers::getNameFromMoabType( const moab::EntityType moab_type )
{
    std::string name;
    switch ( moab_type )
    {
    case moab::MBVERTEX:
        name = "Vertex";
        break;
    case moab::MBEDGE:
        name = "Edge";
        break;
    case moab::MBTRI:
        name = "Triangle";
        break;
    case moab::MBQUAD:
        name = "Quadrilateral";
        break;
    case moab::MBPOLYGON:
        name = "Polygon";
        break;
    case moab::MBTET:
        name = "Tetrahedron";
        break;
    case moab::MBPYRAMID:
        name = "Pyramid";
        break;
    case moab::MBPRISM:
        name = "Prism";
        break;
    case moab::MBKNIFE:
        name = "Knife";
        break;
    case moab::MBHEX:
        name = "Hexahedron";
        break;
    case moab::MBPOLYHEDRON:
        name = "Polyhedron";
        break;
    case moab::MBENTITYSET:
        name = "Entity Set";
        break;
    default:
        name = "Not an entity type";
        break;
    }
    return name;
}

//---------------------------------------------------------------------------//
// Get the coordinates of the entity nodes in canonical order.
void MoabHelpers::getEntityNodeCoordinates(
    const moab::EntityHandle &moab_entity,
    const Teuchos::Ptr<moab::ParallelComm> &moab_mesh,
    Teuchos::Array<double> &coordinates )
{
    const moab::EntityHandle *entity_nodes;
    int num_nodes = 0;
    std::vector<moab::EntityHandle> storage;
    DTK_CHECK_ERROR_CODE( moab_mesh->get_moab()->get_connectivity(
        moab_entity, entity_nodes, num_nodes, false, &storage ) );

    coordinates.resize( 3 * num_nodes );
    DTK_CHECK_ERROR_CODE( moab_mesh->get_moab()->get_coords(
        entity_nodes, num_nodes, coordinates.getRawPtr() ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabHelpers.cpp
//---------------------------------------------------------------------------//
