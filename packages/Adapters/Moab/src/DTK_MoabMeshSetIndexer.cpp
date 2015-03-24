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
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MoabMeshSetIndexer::MoabMeshSetIndexer( 
    const Teuchos::RCP<moab::ParallelComm>& moab_mesh )
{
    moab::EntityHandle root_set = moab_mesh->get_moab()->get_root_set();
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
moab::EntityHandle MoabMeshSetIndexer::getMeshSetFromIndex( const int index ) const
{
    DTK_REQUIRE( d_index_to_handle_map.count(index) );
    return d_index_to_handle_map.find( index )->second;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabMeshSetIndexer.cpp
//---------------------------------------------------------------------------//
