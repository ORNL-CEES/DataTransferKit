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
 * \file DTK_ElementTree_impl.hpp
 * \author Stuart R. Slattery
 * \brief ElementTree definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ELEMENTTREE_IMPL_HPP
#define DTK_ELEMENTTREE_IMPL_HPP

#include <cmath>

#include "DTK_DBC.hpp"
#include "DTK_CellTopologyFactory.hpp"
#include "DTK_TopologyTools.hpp"

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param mesh The mesh over which to build the kD-tree. 
 */
template<class Mesh>
ElementTree<Mesh>::ElementTree( const Teuchos::RCP<MeshManager<Mesh> >& mesh )
  : d_mesh( mesh )
{
    // Setup the centroid array. These will be interleaved.
    d_element_centroids.resize( d_mesh->dim() * d_mesh->localNumElements() );

    // Add the centroids from each block.
    int num_blocks = d_mesh->getNumBlocks();
    d_dcells.resize( num_blocks );
    Teuchos::RCP<Mesh> block;
    int block_num_elements = 0;
    int block_num_vertices = 0;
    int block_verts_per_elem = 0;
    Teuchos::ArrayRCP<const double> block_vertex_coords;
    Teuchos::ArrayRCP<const GlobalOrdinal> block_connectivity;
    Teuchos::RCP<shards::CellTopology> block_topology;
    GlobalOrdinal vertex_gid = 0;
    int vertex_lid = 0;
    int elem_lid = 0;
    Teuchos::Array<int> centroid_array_dims(3);
    centroid_array_dims[0] = 1;
    centroid_array_dims[1] = 1;
    centroid_array_dims[2] = mesh->dim();
    for ( int b = 0; b < num_blocks; ++b )
    {
	// Get the current block.
	block = d_mesh->getBlock( b );

	// Make a discretization cell for this block.
	block_verts_per_elem = MT::verticesPerElement( *block );
	block_topology = CellTopologyFactory::create(
	    MT::elementTopology(*block), block_verts_per_elem );
	d_dcells[b] = IntrepidCell<Intrepid::FieldContainer<double> >(
	    *block_topology, 1 );

	// Get the reference cell center for this topology.
	Intrepid::FieldContainer<double> 
	    ref_center( 1, d_dcells[b].getSpatialDimension() );
	TopologyTools::referenceCellCenter( *block_topology, ref_center );

	// Get the centroid of each element in the block.
	block_num_elements = MeshTools<Mesh>::numElements( *block );
	block_vertex_coords = MeshTools<Mesh>::coordsView( *block );
	block_connectivity = MeshTools<Mesh>::connectivityView( *block );
	block_num_vertices = MeshTools<Mesh>::numVertices( *block );
	Intrepid::FieldContainer<double> node_coords( 
	    1, block_verts_per_elem, d_mesh->dim() );
	for ( int e = 0; e < block_num_elements; ++e, ++elem_lid )
	{
	    // Extract the element node coordinates.
	    for ( int n = 0; n < block_verts_per_elem; ++n )
	    {
		vertex_gid = block_connectivity[ n*block_num_elements + e ];
		vertex_lid = d_mesh->getLocalVertexId( b, vertex_gid );

		for ( int d = 0; d < d_mesh->dim(); ++d )
		{
		    node_coords( 0, n, d ) = 
			block_vertex_coords[ block_num_vertices*d + vertex_lid ];
		}
	    }

	    // Set the state of the discretization cell with the coordinates.
	    d_dcells[b].setCellNodeCoordinates( node_coords );

	    // Map the reference cell center to the physical frame of this
	    // element to construct the centroid.
	    Intrepid::FieldContainer<double> centroid( 
		centroid_array_dims, 
		&d_element_centroids[ elem_lid*d_mesh->dim() ] );
	    d_dcells[b].mapToCellPhysicalFrame( ref_center, centroid );
	}
    }

    // Build a cloud search.
    unsigned leaf_size = std::log( d_mesh->localNumElements() ) + 1;
    leaf_size *= leaf_size;
    d_tree = createSearchTree( d_mesh->dim(), d_element_centroids(), leaf_size );
    DTK_ENSURE( Teuchos::nonnull(d_tree) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename GlobalOrdinal>
ElementTree<GlobalOrdinal>::~ElementTree()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in the tree. Return false if we didn't find it in the
 * tree.
 *
 * \param coords Point coordinates to locate in the tree. Point dimensions
 * less than or equal to 3 are valid but the point most be the same dimension
 * as the tree.
 *
 * \param element The global ordinal of the client element the point was found
 * in. This global ordinal is not valid if this function returns false.
 *
 * \param tolerance Absolute tolerance for point searching. Will be used when
 * checking the reference cell ( and is therefore absolute ).
 *
 * \return Return true if the point was found in the kD-tree, false if not.
 */
template<class Mesh>
bool ElementTree<Mesh>::findPoint( 
    const Teuchos::ArrayView<const double>& coords,
    GlobalOrdinal& element,
    double tolerance )
{
    // Wrap the point in a field container.
    Teuchos::Array<int> point_array_dims(2);
    point_array_dims[0] = 1;
    point_array_dims[1] = d_mesh->dim();
    Intrepid::FieldContainer<double> point( 
	point_array_dims, const_cast<double*>(coords.getRawPtr()) );

    // Find the leaf of nearest neighbors.
    unsigned leaf_size = std::log( d_mesh->localNumElements() ) + 1;
    unsigned num_neighbors = leaf_size*leaf_size;
    Teuchos::Array<unsigned> neighbors = 
	d_tree->nnSearch( coords, num_neighbors );

    // Check the leaf for point inclusion.
    int block_id = 0;
    int elem_id = 0;
    Teuchos::RCP<Mesh> block;
    Teuchos::ArrayRCP<const double> block_vertex_coords;
    Teuchos::ArrayRCP<const GlobalOrdinal> block_element_gids;
    Teuchos::ArrayRCP<const GlobalOrdinal> block_connectivity;
    int block_verts_per_elem = 0;
    int block_num_elements = 0;
    int block_num_vertices = 0;
    GlobalOrdinal vertex_gid = 0;
    int vertex_lid = 0;
    bool point_in_element = false;

    for ( unsigned n = 0; n < neighbors.size(); ++n )
    {
	// Get the block and element id for this neighbor resides in.
	d_mesh->getLocalElementIds( neighbors[n], block_id, elem_id );
	block = d_mesh->getBlock(block_id);

	// Get a view of the connectivity and vertex coordinates of that
	// block.
	block_vertex_coords = MeshTools<Mesh>::coordsView( *block );
	block_connectivity = MeshTools<Mesh>::connectivityView( *block );
	block_verts_per_elem = MT::verticesPerElement( *block );
	block_num_elements = MeshTools<Mesh>::numElements( *block );
	block_num_vertices = MeshTools<Mesh>::numVertices( *block );

	// Extract the element node coordinates.
	Intrepid::FieldContainer<double> node_coords( 
	    1, block_verts_per_elem, d_mesh->dim() );
	for ( int i = 0; i < block_verts_per_elem; ++i )
	{
	    vertex_gid = block_connectivity[ i*block_num_elements + elem_id ];
	    vertex_lid = d_mesh->getLocalVertexId( block_id, vertex_gid );

	    for ( int d = 0; d < d_mesh->dim(); ++d )
	    {
		node_coords( 0, i, d ) = 
		    block_vertex_coords[ block_num_vertices*d + vertex_lid ];
	    }
	}

	// Set the state of the block discretization cell with the
	// coordinates.
	d_dcells[block_id].setCellNodeCoordinates( node_coords );

	// Check for point inclusion.
	point_in_element = 
	    d_dcells[block_id].pointInPhysicalCell( point, tolerance );

	// If we found the point in one of the elements set the element id to
	// this element and return true.
	if ( point_in_element )
	{
	    block_element_gids = MeshTools<Mesh>::elementsView(*block);
	    element = block_element_gids[elem_id];
	    return true;
	}
    }

    // If we got here we did not find the point in any element.
    element = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
    return false;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ELEMENTTREE_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ElementTree_impl.hpp
//---------------------------------------------------------------------------//

