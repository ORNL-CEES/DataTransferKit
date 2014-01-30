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
 * \file DTK_MeshManager_def.hpp
 * \author Stuart R. Slattery
 * \brief Mesh manager definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHMANAGER_DEF_HPP
#define DTK_MESHMANAGER_DEF_HPP

#include <algorithm>
#include <limits>

#include "DTK_MeshTypes.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_DBC.hpp"
#include "DataTransferKit_config.hpp"

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
template<class Mesh>
MeshManager<Mesh>::MeshManager( 
    const Teuchos::ArrayRCP<RCP_Mesh>& mesh_blocks,
    const RCP_Comm& comm, const int dim )
    : d_mesh_blocks( mesh_blocks )
    , d_comm( comm )
    , d_dim( dim )
    , d_active_vertices( d_mesh_blocks.size() )
    , d_active_elements( d_mesh_blocks.size() )
{
    // If we're checking with Design-by-Contract, validate the mesh to the
    // domain model.
#if HAVE_DTK_DBC
    validate();
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh>
MeshManager<Mesh>::~MeshManager()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local number of elements in the mesh across all mesh blocks.
 *
 * \return The local number of elements across all mesh blocks.
 */
template<class Mesh>
typename MeshManager<Mesh>::GlobalOrdinal 
MeshManager<Mesh>::localNumElements() const
{
    GlobalOrdinal local_num_elements = 0;
    BlockIterator block_iterator;
    for ( block_iterator = d_mesh_blocks.begin();
	  block_iterator != d_mesh_blocks.end();
	  ++block_iterator )
    {
	local_num_elements += 
	    MeshTools<Mesh>::numElements( *(*block_iterator) );
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
template<class Mesh>
typename MeshManager<Mesh>::GlobalOrdinal 
MeshManager<Mesh>::globalNumElements() const
{
    GlobalOrdinal local_num_elements = localNumElements();
    GlobalOrdinal global_num_elements = 0;
    Teuchos::reduceAll<int,GlobalOrdinal>( *d_comm,
					   Teuchos::REDUCE_SUM,
					   1,
					   &local_num_elements,
					   &global_num_elements );
    return global_num_elements;
}

//---------------------------------------------------------------------------//
/*!    
 * \brief Compute the global bounding box around the entire mesh (all
 * blocks). 
 *
 * \return The global bounding box over all mesh blocks.
 */
template<class Mesh>
BoundingBox MeshManager<Mesh>::globalBoundingBox()
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
    BlockIterator block_iterator;
    for ( block_iterator = d_mesh_blocks.begin();
	  block_iterator != d_mesh_blocks.end();
	  ++block_iterator )
    {
	// If the mesh block is empty, do nothing.
	if ( MeshTools<Mesh>::numVertices( *(*block_iterator) ) > 0 )
	{
	    block_box =	MeshTools<Mesh>::globalBoundingBox( 
		*(*block_iterator), d_comm );

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
 * \brief Validate the mesh to the domain model.
 */
template<class Mesh>
void MeshManager<Mesh>::validate()
{
    // Check that the mesh is of a valid dimension.
    DTK_REQUIRE( 0 <= d_dim && d_dim <= 3 );

    // Check that the mesh dimension is the same on every node.
    Teuchos::Array<int> local_dims( d_comm->getSize(), 0 );
    Teuchos::Array<int> local_dims_copy( d_comm->getSize(), 0 );
    local_dims[ d_comm->getRank() ] = d_dim;
    Teuchos::reduceAll<int,int>( *d_comm, Teuchos::REDUCE_SUM,
				 local_dims.size(),
				 &local_dims[0], &local_dims_copy[0] ); 
    Teuchos::Array<int>::iterator unique_bound;
    std::sort( local_dims_copy.begin(), local_dims_copy.end() );
    unique_bound = std::unique( local_dims_copy.begin(), local_dims_copy.end() );
    int unique_dim = std::distance( local_dims_copy.begin(), unique_bound );
    DTK_REQUIRE( 1 == unique_dim );
    local_dims_copy.clear();

    // Check that the same number of blocks have been defined on every node.
    Teuchos::Array<int> local_blocks( d_comm->getSize(), 0 );
    Teuchos::Array<int> local_blocks_copy( d_comm->getSize(), 0 );
    local_blocks[ d_comm->getRank() ] = getNumBlocks();
    Teuchos::reduceAll<int,int>( *d_comm, Teuchos::REDUCE_SUM,
				 local_blocks.size(),
				 &local_blocks[0], &local_blocks_copy[0] ); 
    std::sort( local_blocks_copy.begin(), local_blocks_copy.end() );
    unique_bound = std::unique( local_blocks_copy.begin(), local_blocks_copy.end() );
    int unique_blocks = std::distance( local_blocks_copy.begin(), unique_bound );
    DTK_REQUIRE( 1 == unique_blocks );
    local_blocks_copy.clear();

    // Check the mesh blocks.
    BlockIterator block_iterator;
    for ( block_iterator = d_mesh_blocks.begin();
	  block_iterator != d_mesh_blocks.end();
	  ++block_iterator )
    {
	// Check that the block vertices are the same dimension as the mesh.
	DTK_REQUIRE( d_dim == MT::vertexDim( *(*block_iterator) ) );

	// Check that the coordinate dimension is the same as the mesh
	// dimension.
	GlobalOrdinal num_vertices = 
	    MeshTools<Mesh>::numVertices( *(*block_iterator) );
	GlobalOrdinal num_coords = std::distance( 
	    MT::coordsBegin( *(*block_iterator) ), 
	    MT::coordsEnd( *(*block_iterator) ) );
	if ( num_vertices > 0 )
	{
	    DTK_REQUIRE( num_coords / num_vertices 
			      == Teuchos::as<GlobalOrdinal>(d_dim) );
	}
	
	// Check that the element topology is valid for the given dimension.
	if ( d_dim == 0 )
	{
	    DTK_REQUIRE( MT::elementTopology( *(*block_iterator) ) 
			      == DTK_VERTEX );
	}
	else if ( d_dim == 1 )
	{
	    DTK_REQUIRE( MT::elementTopology( *(*block_iterator) ) == 
			      DTK_LINE_SEGMENT );
	}
	else if ( d_dim == 2 )
	{
	    DTK_REQUIRE( MT::elementTopology( *(*block_iterator) ) == 
			      DTK_TRIANGLE ||
			      MT::elementTopology( *(*block_iterator) ) == 
			      DTK_QUADRILATERAL );
	}
	else if ( d_dim == 3 )
	{
	    DTK_REQUIRE( MT::elementTopology( *(*block_iterator) ) == 
			      DTK_TETRAHEDRON ||
			      MT::elementTopology( *(*block_iterator) ) == 
			      DTK_HEXAHEDRON ||
			      MT::elementTopology( *(*block_iterator) ) == 
			      DTK_PYRAMID ||
			      MT::elementTopology( *(*block_iterator) ) == 
			      DTK_WEDGE );
	}

	// Check that this block has the same topology on all nodes.
	Teuchos::Array<int> local_topo( d_comm->getSize(), 0 );
	Teuchos::Array<int> local_topo_copy( d_comm->getSize(), 0 );
	local_topo[ d_comm->getRank() ] = 
	    Teuchos::as<int>(MT::elementTopology( *(*block_iterator) ));
	Teuchos::reduceAll<int,int>( *d_comm, Teuchos::REDUCE_SUM,
				     local_topo.size(),
				     &local_topo[0], &local_topo_copy[0] ); 
	std::sort( local_topo_copy.begin(), local_topo_copy.end() );
	unique_bound = std::unique( local_topo_copy.begin(), local_topo_copy.end() );
	int unique_topo = std::distance( local_topo_copy.begin(), unique_bound );
	DTK_REQUIRE( 1 == unique_topo );
	local_topo_copy.clear();

	// Check that the element handles are of a value less than the numeric
	// limit of the ordinal type. This is an invalid element handle.
	typename MT::const_element_iterator element_iterator;
	for ( element_iterator = MT::elementsBegin( *(*block_iterator) );
	      element_iterator != MT::elementsEnd( *(*block_iterator) );
	      ++element_iterator )
	{
	    DTK_REQUIRE( *element_iterator < 
			      std::numeric_limits<GlobalOrdinal>::max() );
	}
	
	// Check that the connectivity size is the same as the number of
	// vertices per element.
	GlobalOrdinal num_elements =
	    MeshTools<Mesh>::numElements( *(*block_iterator) );
	GlobalOrdinal num_conn = std::distance( 
	    MT::connectivityBegin( *(*block_iterator) ),
	    MT::connectivityEnd( *(*block_iterator) ) );
	if ( num_elements > Teuchos::as<GlobalOrdinal>(0) )
	{
	    DTK_REQUIRE( num_conn / num_elements ==
			      Teuchos::as<GlobalOrdinal>(
				  MT::verticesPerElement(*(*block_iterator))) );
	}

	// Check that the size of the permutation vector is the same as the
	// number of vertices per element.
	int num_permutation = std::distance(
	    MT::permutationBegin( *(*block_iterator) ),
	    MT::permutationEnd( *(*block_iterator) ) );
	DTK_REQUIRE( MT::verticesPerElement( *(*block_iterator) ) ==
			  num_permutation );

	// Check that the permutation vector contains unique values.
	Teuchos::Array<int> permutation( num_permutation );
	std::copy( MT::permutationBegin( *(*block_iterator) ),
		   MT::permutationEnd( *(*block_iterator) ), 
		   permutation.begin() );
	Teuchos::Array<int>::iterator permutation_bound;
	std::sort( permutation.begin(), permutation.end() );
	permutation_bound = std::unique( permutation.begin(),
					 permutation.end() );
	int unique_permutation = std::distance( permutation.begin(),
						permutation_bound );
	DTK_REQUIRE( MT::verticesPerElement( *(*block_iterator) ) ==
			  unique_permutation );

	// Check that the permutation vector contains value less than its
	// size. This implies that we shouldn't get more vertices in the
	// connectivity list than are needed to build the linear element.
	typename MT::const_permutation_iterator permutation_it;
	for ( permutation_it = MT::permutationBegin( *(*block_iterator) );
	      permutation_it != MT::permutationEnd( *(*block_iterator) );
	      ++permutation_it )
	{
	    DTK_REQUIRE( *permutation_it < 
			      MT::verticesPerElement( *(*block_iterator) ) );
	}

	d_comm->barrier();
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESHMANAGER_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshManager_def.hpp
//---------------------------------------------------------------------------//


