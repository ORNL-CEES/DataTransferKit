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
 * \brief DTK_MoabMeshSetIndexer.cpp
 * \author Stuart R. Slattery
 * \brief Moab mesh set indexer.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "DTK_MoabMeshSetIndexer.hpp"
#include "DTK_MoabHelpers.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MoabMeshSetIndexer::MoabMeshSetIndexer( 
    const Teuchos::RCP<moab::ParallelComm>& moab_mesh )
    : d_gid_map( 4 )
{
    // DTK needs unique global ids from MOAB. Assign global ids to the root
    // set if needed.
    int dimension = 0;
    DTK_CHECK_ERROR_CODE(
	moab_mesh->get_moab()->get_dimension( dimension )
	);
    moab::EntityHandle root_set = moab_mesh->get_moab()->get_root_set();
    DTK_CHECK_ERROR_CODE(
	moab_mesh->assign_global_ids( root_set,
				      dimension,
				      1,
				      false,
				      true,
				      false )
	);

    // Index the entity sets.
    int root_index = 0;
    d_handle_to_index_map.insert( std::make_pair(root_set, root_index) );
    d_index_to_handle_map.insert( std::make_pair(root_index, root_set) );

    std::vector<moab::EntityHandle> mesh_sets;
    DTK_CHECK_ERROR_CODE(
	moab_mesh->get_moab()->get_contained_meshsets( root_set, mesh_sets, 0 )
	);

    int num_sets = mesh_sets.size();
    for ( int i = 0; i < num_sets; ++i )
    {
	d_handle_to_index_map.insert( std::make_pair(mesh_sets[i], i+1) );
	d_index_to_handle_map.insert( std::make_pair(i+1, mesh_sets[i]) );
    }

    // Map global ids to entities.
    int space_dim = 0;
    DTK_CHECK_ERROR_CODE(
	moab_mesh->get_moab()->get_dimension( space_dim )
	);
    std::vector<moab::EntityHandle> dim_entities;
    std::vector<moab::EntityHandle>::const_iterator entity_it;
    std::vector<EntityId> gid_data;
    std::vector<EntityId>::const_iterator gid_it;
    for ( int d = 0; d < space_dim+1; ++d )
    {
	// Get the dimension entities.
	dim_entities.clear();
	DTK_CHECK_ERROR_CODE(
	    moab_mesh->get_moab()->get_entities_by_dimension(
		0, d, dim_entities )
	    );

	// Get the ids.
	gid_data.resize( dim_entities.size() );
	MoabHelpers::getGlobalIds( *moab_mesh,
				   dim_entities.data(),
				   dim_entities.size(),
				   gid_data.data() );

	// Map the gids to the enties.
	for ( entity_it = dim_entities.begin(),
		 gid_it = gid_data.begin();
	      entity_it != dim_entities.end();
	      ++entity_it, ++gid_it )
	{
	    d_gid_map[d].emplace( *gid_it, *entity_it );
	}
    }
}

//---------------------------------------------------------------------------//
// Given an entity set handle, get the integer index in the mesh.
int MoabMeshSetIndexer::getIndexFromMeshSet( 
    const moab::EntityHandle mesh_set ) const
{
    DTK_REQUIRE( d_handle_to_index_map.count(mesh_set) );
    return d_handle_to_index_map.find( mesh_set )->second;
}

//---------------------------------------------------------------------------//
// Given an integer index, get the entity set handle.
moab::EntityHandle
MoabMeshSetIndexer::getMeshSetFromIndex( const int index ) const
{
    DTK_REQUIRE( d_index_to_handle_map.count(index) );
    return d_index_to_handle_map.find( index )->second;
}

//---------------------------------------------------------------------------//
// Given a global id and topological, get its entity.
moab::EntityHandle MoabMeshSetIndexer::getEntityFromGlobalId(
    const EntityId id,
    const int topological_dimension ) const
{
    DTK_REQUIRE( d_gid_map[topological_dimension].count(id) );
    return d_gid_map[topological_dimension].find( id )->second;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabMeshSetIndexer.cpp
//---------------------------------------------------------------------------//
