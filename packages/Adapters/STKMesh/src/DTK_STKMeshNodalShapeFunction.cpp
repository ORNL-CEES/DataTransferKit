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
 * \brief DTK_STKMeshNodalShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief Nodal shape function implementation for STK mesh.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshNodalShapeFunction.hpp"
#include "DTK_DBC.hpp"
#include "DTK_STKMeshHelpers.hpp"

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
STKMeshNodalShapeFunction::STKMeshNodalShapeFunction(
    const Teuchos::RCP<stk::mesh::BulkData> &bulk_data )
    : d_bulk_data( bulk_data )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Given an entity, get the ids of the degrees of freedom in the vector space
// supporting its shape function.
void STKMeshNodalShapeFunction::entitySupportIds(
    const Entity &entity, Teuchos::Array<SupportId> &support_ids ) const
{
    // Extract the stk entity.
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );
    const stk::mesh::EntityRank entity_rank =
        d_bulk_data->entity_rank( stk_entity );

    // If the entity is a node, return the id of the node as the support id.
    if ( stk::topology::NODE_RANK == entity_rank )
    {
        support_ids.assign( 1, d_bulk_data->identifier( stk_entity ) );
    }

    // Otherwise get the ids of the nodes supporting the entity.
    else
    {
        const stk::mesh::Entity *begin = d_bulk_data->begin_nodes( stk_entity );
        const stk::mesh::Entity *end = d_bulk_data->end_nodes( stk_entity );

        // Extract the node ids as the support ids.
        int num_nodes = std::distance( begin, end );
        support_ids.resize( num_nodes );
        for ( int n = 0; n < num_nodes; ++n )
        {
            support_ids[n] = d_bulk_data->identifier( begin[n] );
        }
    }
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void STKMeshNodalShapeFunction::evaluateValue(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<double> &values ) const
{
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );

    shards::CellTopology entity_topo = stk::mesh::get_cell_topology(
        d_bulk_data->bucket( stk_entity ).topology() );

    d_intrepid_shape.evaluateValue( entity_topo, reference_point, values );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void STKMeshNodalShapeFunction::evaluateGradient(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<Teuchos::Array<double>> &gradients ) const
{
    const stk::mesh::Entity &stk_entity =
        STKMeshHelpers::extractEntity( entity );

    shards::CellTopology entity_topo = stk::mesh::get_cell_topology(
        d_bulk_data->bucket( stk_entity ).topology() );

    d_intrepid_shape.evaluateGradient( entity_topo, reference_point,
                                       gradients );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshNodalShapeFunction.cpp
//---------------------------------------------------------------------------//
