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
 * \brief DTK_STKMeshHelpers_impl.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh helpers.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHHELPERS_IMPL_HPP
#define DTK_STKMESHHELPERS_IMPL_HPP

#include <limits>
#include <unordered_set>

#include "DTK_DBC.hpp"

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Extract the coordinates of multiple entities and return an array ordered as
// (C,N,D).
template <class FieldType>
Intrepid::FieldContainer<double> STKMeshHelpers::extractEntityNodeCoordinates(
    const Teuchos::Array<stk::mesh::Entity> &stk_entities,
    const stk::mesh::BulkData &bulk_data, const int space_dim )
{
    // Cast the field.
    const stk::mesh::FieldBase *coord_field_base =
        bulk_data.mesh_meta_data().coordinate_field();
    const stk::mesh::Field<double, FieldType> *coord_field =
        dynamic_cast<const stk::mesh::Field<double, FieldType> *>(
            coord_field_base );

    // Allocate the coordinate array.
    int num_cells = stk_entities.size();
    int num_nodes = 0;
    stk::mesh::EntityRank stk_rank = stk::topology::INVALID_RANK;
    if ( num_cells > 0 )
    {
        stk_rank = bulk_data.entity_rank( stk_entities[0] );
        if ( stk::topology::NODE_RANK == stk_rank )
        {
            num_nodes = 1;
        }
        else
        {
            const stk::mesh::Entity *begin =
                bulk_data.begin_nodes( stk_entities[0] );
            const stk::mesh::Entity *end =
                bulk_data.end_nodes( stk_entities[0] );
            num_nodes = std::distance( begin, end );
        }
    }
    Intrepid::FieldContainer<double> coords( num_cells, num_nodes, space_dim );

    // Extract the coordinates.
    double *node_coords = 0;
    for ( int c = 0; c < num_cells; ++c )
    {
        if ( stk::topology::NODE_RANK == stk_rank )
        {
            node_coords =
                stk::mesh::field_data( *coord_field, stk_entities[c] );
            for ( int d = 0; d < space_dim; ++d )
            {
                coords( c, 0, d ) = node_coords[d];
            }
        }
        else
        {
            const stk::mesh::Entity *begin =
                bulk_data.begin_nodes( stk_entities[c] );
            DTK_REMEMBER( const stk::mesh::Entity *end =
                              bulk_data.end_nodes( stk_entities[c] ) );
            DTK_CHECK( std::distance( begin, end ) == num_nodes );
            for ( int n = 0; n < num_nodes; ++n )
            {
                node_coords = stk::mesh::field_data( *coord_field, begin[n] );
                for ( int d = 0; d < space_dim; ++d )
                {
                    coords( c, n, d ) = node_coords[d];
                }
            }
        }
    }

    return coords;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_STKMESHHELPERS_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshHelpers_impl.hpp
//---------------------------------------------------------------------------//
