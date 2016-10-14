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
 * \brief DTK_STKMeshEntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshEntityLocalMap.hpp"
#include "DTK_DBC.hpp"
#include "DTK_IntrepidCellLocalMap.hpp"
#include "DTK_ProjectionPrimitives.hpp"
#include "DTK_STKMeshHelpers.hpp"

#include <Intrepid_FieldContainer.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
STKMeshEntityLocalMap::STKMeshEntityLocalMap(
    const Teuchos::RCP<stk::mesh::BulkData> &bulk_data )
    : d_bulk_data( bulk_data )
    , d_inclusion_tol( 1.0e-6 )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Set parameters for mapping.
void STKMeshEntityLocalMap::setParameters(
    const Teuchos::ParameterList &parameters )
{
    if ( parameters.isParameter( "Point Inclusion Tolerance" ) )
    {
        d_inclusion_tol = parameters.get<double>( "Point Inclusion Tolerance" );
    }
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
double STKMeshEntityLocalMap::measure( const Entity &entity ) const
{
    // Get the STK entity and its topology.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    stk::mesh::EntityRank rank = d_bulk_data->entity_rank( stk_entity );
    shards::CellTopology entity_topo =
        STKMeshHelpers::getShardsTopology( stk_entity, *d_bulk_data );

    // Get the STK entity coordinates.
    Intrepid::FieldContainer<double> entity_coords =
        STKMeshHelpers::getEntityNodeCoordinates(
            Teuchos::Array<stk::mesh::Entity>( 1, stk_entity ), *d_bulk_data );

    // Compute the measure of the element.
    if ( rank == stk::topology::ELEM_RANK )
    {
        return IntrepidCellLocalMap::measure( entity_topo, entity_coords );
    }

    // Compute the measure of the face.
    else if ( rank == stk::topology::FACE_RANK )
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
        return -1.0;
    }

    // Check for unsupported ranks.
    else
    {
        bool bad_rank = true;
        DTK_INSIST( !bad_rank );
        return -1.0;
    }
    return -1.0;
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void STKMeshEntityLocalMap::centroid(
    const Entity &entity, const Teuchos::ArrayView<double> &centroid ) const
{
    // Get the STK entity.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    stk::mesh::EntityRank rank = d_bulk_data->entity_rank( stk_entity );

    // Extract the centroid of the element.
    if ( rank == stk::topology::ELEM_RANK )
    {
        shards::CellTopology entity_topo =
            STKMeshHelpers::getShardsTopology( stk_entity, *d_bulk_data );
        Intrepid::FieldContainer<double> entity_coords =
            STKMeshHelpers::getEntityNodeCoordinates(
                Teuchos::Array<stk::mesh::Entity>( 1, stk_entity ),
                *d_bulk_data );
        IntrepidCellLocalMap::centroid( entity_topo, entity_coords, centroid );
    }

    // Extract the centroid of the face.
    else if ( rank == stk::topology::FACE_RANK )
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
    }

    // The centroid of a node is the node coordinates.
    else if ( rank == stk::topology::NODE_RANK )
    {
        Intrepid::FieldContainer<double> entity_coords =
            STKMeshHelpers::getEntityNodeCoordinates(
                Teuchos::Array<stk::mesh::Entity>( 1, stk_entity ),
                *d_bulk_data );
        centroid.assign( entity_coords.getData()() );
    }

    // Check for unsupported ranks.
    else
    {
        bool bad_rank = true;
        DTK_INSIST( !bad_rank );
    }
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance.
bool STKMeshEntityLocalMap::isSafeToMapToReferenceFrame(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &physical_point ) const
{
    // Get the STK entity.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    stk::mesh::EntityRank rank = d_bulk_data->entity_rank( stk_entity );

    // If we have an element, use the default implementation.
    if ( rank == stk::topology::ELEM_RANK )
    {
        return EntityLocalMap::isSafeToMapToReferenceFrame( entity,
                                                            physical_point );
    }

    // If we have a face, perform the projection safeguard.
    else if ( rank == stk::topology::FACE_RANK )
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
        return false;
    }

    // Check for unsupported ranks.
    else
    {
        bool bad_rank = true;
        DTK_INSIST( !bad_rank );
        return false;
    }

    // Return true to indicate successful mapping. Catching Intrepid errors
    // and returning false is a possibility here.
    return true;
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
bool STKMeshEntityLocalMap::mapToReferenceFrame(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &physical_point,
    const Teuchos::ArrayView<double> &reference_point ) const
{
    // Get the STK entity.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    stk::mesh::EntityRank rank = d_bulk_data->entity_rank( stk_entity );

    // Use the cell to perform the element mapping.
    if ( rank == stk::topology::ELEM_RANK )
    {
        shards::CellTopology entity_topo =
            STKMeshHelpers::getShardsTopology( stk_entity, *d_bulk_data );
        Intrepid::FieldContainer<double> entity_coords =
            STKMeshHelpers::getEntityNodeCoordinates(
                Teuchos::Array<stk::mesh::Entity>( 1, stk_entity ),
                *d_bulk_data );
        IntrepidCellLocalMap::mapToReferenceFrame(
            entity_topo, entity_coords, physical_point, reference_point );
    }

    // Use the side cell to perform the face mapping.
    else if ( rank == stk::topology::FACE_RANK )
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
        return false;
    }

    // Check for unsupported ranks.
    else
    {
        bool bad_rank = true;
        DTK_INSIST( !bad_rank );
    }

    // Return true to indicate successful mapping. Catching Intrepid errors
    // and returning false is a possibility here.
    return true;
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool STKMeshEntityLocalMap::checkPointInclusion(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    // Get the STK entity and its topology.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    stk::mesh::EntityRank rank = d_bulk_data->entity_rank( stk_entity );
    shards::CellTopology entity_topo =
        STKMeshHelpers::getShardsTopology( stk_entity, *d_bulk_data );

    // Check point inclusion in the element.
    if ( rank == stk::topology::ELEM_RANK )
    {
        return IntrepidCellLocalMap::checkPointInclusion(
            entity_topo, reference_point, d_inclusion_tol );
    }

    // Check point inclusion in the face.
    else if ( rank == stk::topology::FACE_RANK )
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
        return false;
    }

    // Check for unsupported ranks.
    else
    {
        bool bad_rank = true;
        DTK_INSIST( !bad_rank );
        return false;
    }

    return false;
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void STKMeshEntityLocalMap::mapToPhysicalFrame(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &physical_point ) const
{
    // Get the STK entity.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    stk::mesh::EntityRank rank = d_bulk_data->entity_rank( stk_entity );

    // Map from the element.
    if ( rank == stk::topology::ELEM_RANK )
    {
        shards::CellTopology entity_topo =
            STKMeshHelpers::getShardsTopology( stk_entity, *d_bulk_data );
        Intrepid::FieldContainer<double> entity_coords =
            STKMeshHelpers::getEntityNodeCoordinates(
                Teuchos::Array<stk::mesh::Entity>( 1, stk_entity ),
                *d_bulk_data );
        IntrepidCellLocalMap::mapToPhysicalFrame(
            entity_topo, entity_coords, reference_point, physical_point );
    }

    // Map from the face.
    else if ( rank == stk::topology::FACE_RANK )
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
    }

    // Check for unsupported ranks.
    else
    {
        bool bad_rank = true;
        DTK_INSIST( !bad_rank );
    }
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference point.
void STKMeshEntityLocalMap::normalAtReferencePoint(
    const Entity &entity, const Entity &parent_entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &normal ) const
{
    // Get the STK entity.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    stk::mesh::EntityRank rank = d_bulk_data->entity_rank( stk_entity );

    // We can only compute normals for faces.
    if ( rank == stk::topology::FACE_RANK )
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
    }

    // Check for unsupported ranks.
    else
    {
        bool bad_rank = true;
        DTK_INSIST( !bad_rank );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityLocalMap.cpp
//---------------------------------------------------------------------------//
