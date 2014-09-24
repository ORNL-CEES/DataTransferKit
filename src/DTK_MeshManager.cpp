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
 * \file DTK_MeshManager.cpp
 * \author Stuart R. Slattery
 * \brief Mesh manager definition.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <limits>

#include "DTK_MeshManager.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor. If Design-By-Contract is enabled, the constructor will
 * validate the mesh description to the domain model. This requires a few
 * global communications.
 *
 * \param mesh_blocks The blocks that construct this mesh. Each block must
 * have MeshTraits.
 *
 * \param comm The communicator the mesh is defined over.
 *
 * \param dim The mesh dimension.
 */
MeshManager::MeshManager( 
    const Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> >& mesh_blocks,
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
    const int dim )
    : d_mesh_blocks( mesh_blocks )
    , d_comm( comm )
    , d_dim( dim )
    , d_active_vertices( d_mesh_blocks.size() )
    , d_active_elements( d_mesh_blocks.size() )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
MeshManager::~MeshManager()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local number of elements in the mesh across all mesh blocks.
 *
 * \return The local number of elements across all mesh blocks.
 */
MeshId MeshManager::localNumElements() const
{
    MeshId local_num_elements = 0;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = d_mesh_blocks.begin();
	  block_iterator != d_mesh_blocks.end();
	  ++block_iterator )
    {
	local_num_elements += (*block_iterator)->elementIds().size();
    }
    return local_num_elements;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global number of elements in the mesh across all mesh
 * blocks. 
 *
 * \return The global number of elements in the mesh across all mesh blocks.
 */
MeshId MeshManager::globalNumElements() const
{
    MeshId local_num_elements = localNumElements();
    MeshId global_num_elements = 0;
    Teuchos::reduceAll<int,MeshId>( *d_comm,
				    Teuchos::REDUCE_SUM,
				    1,
				    &local_num_elements,
				    &global_num_elements );
    return global_num_elements;
}

//---------------------------------------------------------------------------//
/*!    
 * \brief Compute the local bounding box around the entire mesh (all
 * blocks). 
 *
 * \return The local bounding box over all mesh blocks.
 */
BoundingBox MeshManager::localBoundingBox()
{
    // Check to see if there are any local elements.
    if ( 0 == localNumElements() )
    {
	return BoundingBox( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    }

    double local_x_min = Teuchos::ScalarTraits<double>::rmax();
    double local_y_min = Teuchos::ScalarTraits<double>::rmax();
    double local_z_min = Teuchos::ScalarTraits<double>::rmax();
    double local_x_max = -Teuchos::ScalarTraits<double>::rmax();
    double local_y_max = -Teuchos::ScalarTraits<double>::rmax();
    double local_z_max = -Teuchos::ScalarTraits<double>::rmax();

    // Get the bounding box for each mesh block.
    Teuchos::Tuple<double,6> box_bounds;
    BoundingBox block_box;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = d_mesh_blocks.begin();
	  block_iterator != d_mesh_blocks.end();
	  ++block_iterator )
    {
	// If the mesh block is empty, do nothing.
	if ( (*block_iterator)->vertexIds().size() > 0 )
	{
	    block_box =	MeshTools::localBoundingBox( *block_iterator );

	    box_bounds = block_box.getBounds();

	    if ( box_bounds[0] < local_x_min )
	    {
		local_x_min = box_bounds[0];
	    }
	    if ( box_bounds[1] < local_y_min )
	    {
		local_y_min = box_bounds[1];
	    }
	    if ( box_bounds[2] < local_z_min )
	    {
		local_z_min = box_bounds[2];
	    }
	    if ( box_bounds[3] > local_x_max )
	    {
		local_x_max = box_bounds[3];
	    }
	    if ( box_bounds[4] > local_y_max )
	    {
		local_y_max = box_bounds[4];
	    }
	    if ( box_bounds[5] > local_z_max )
	    {
		local_z_max = box_bounds[5];
	    }
	}
    }

    return BoundingBox( local_x_min, local_y_min, local_z_min,
			local_x_max, local_y_max, local_z_max );
}

//---------------------------------------------------------------------------//
/*!    
 * \brief Compute the global bounding box around the entire mesh (all
 * blocks). 
 *
 * \return The global bounding box over all mesh blocks.
 */
BoundingBox MeshManager::globalBoundingBox()
{
    double global_x_min = Teuchos::ScalarTraits<double>::rmax();
    double global_y_min = Teuchos::ScalarTraits<double>::rmax();
    double global_z_min = Teuchos::ScalarTraits<double>::rmax();
    double global_x_max = -Teuchos::ScalarTraits<double>::rmax();
    double global_y_max = -Teuchos::ScalarTraits<double>::rmax();
    double global_z_max = -Teuchos::ScalarTraits<double>::rmax();

    // Get the bounding box for each mesh block.
    Teuchos::Tuple<double,6> box_bounds;
    BoundingBox block_box;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = d_mesh_blocks.begin();
	  block_iterator != d_mesh_blocks.end();
	  ++block_iterator )
    {
	// If the mesh block is empty, do nothing.
	if ( (*block_iterator)->vertexIds().size() > 0 )
	{
	    block_box =	MeshTools::globalBoundingBox( 
		*block_iterator, d_comm );

	    box_bounds = block_box.getBounds();

	    if ( box_bounds[0] < global_x_min )
	    {
		global_x_min = box_bounds[0];
	    }
	    if ( box_bounds[1] < global_y_min )
	    {
		global_y_min = box_bounds[1];
	    }
	    if ( box_bounds[2] < global_z_min )
	    {
		global_z_min = box_bounds[2];
	    }
	    if ( box_bounds[3] > global_x_max )
	    {
		global_x_max = box_bounds[3];
	    }
	    if ( box_bounds[4] > global_y_max )
	    {
		global_y_max = box_bounds[4];
	    }
	    if ( box_bounds[5] > global_z_max )
	    {
		global_z_max = box_bounds[5];
	    }
	}
    }

    return BoundingBox( global_x_min, global_y_min, global_z_min,
			global_x_max, global_y_max, global_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given the local id of an element in the mesh, get its block id and
 * its local id in the block.
 */
void MeshManager::buildIndexing()
{
    int num_blocks = d_mesh_blocks.size();
    d_cumulative_elements.resize( num_blocks );
    d_vertex_g2l.resize( num_blocks );
    Teuchos::ArrayRCP<const MeshId> block_vertex_ids;
    typename Teuchos::ArrayRCP<const MeshId>::const_iterator vertex_id_it;
    int vertex_lid = 0;
    MeshId num_elements = 0;
    for ( int b = 0; b < num_blocks; ++b )
    {
	// Element accumulation.
	num_elements = d_mesh_blocks[b]->elementIds().size();
	d_cumulative_elements[b] = 
	    ( b > 0 ) 
	    ? num_elements + d_cumulative_elements[b-1] 
	    : num_elements;

	// Global-to-local vertex id mapping.
	block_vertex_ids = d_mesh_blocks[b]->vertexIds();
	vertex_lid = 0;
	for ( vertex_id_it = block_vertex_ids.begin();
	      vertex_id_it != block_vertex_ids.end();
	      ++vertex_id_it, ++vertex_lid )
	{
	    d_vertex_g2l[b][*vertex_id_it] = vertex_lid;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given the local id of an element in the mesh, get its block id and
 * its local id in the block.
 */
void MeshManager::getLocalElementIds( const int local_elem_id,
					    int& block_id,
					    int& block_elem_id ) const
{
    block_id = std::distance(
	d_cumulative_elements.begin(),
	std::upper_bound( d_cumulative_elements.begin(),
			  d_cumulative_elements.end(),
			  local_elem_id) );

    block_elem_id = (block_id > 0)
		    ? local_elem_id - d_cumulative_elements[block_id-1]
		    : local_elem_id;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given the block id of a vertex and its global id, get its local id
 * in that block
 */
int MeshManager::getLocalVertexId( 
    const int block_id, const int vertex_gid ) const
{
    return d_vertex_g2l[block_id].find( vertex_gid )->second;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MeshManager.cpp
//---------------------------------------------------------------------------//


