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

#include <vector>
#include <cmath>

#include "DTK_TopologyTools.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_as.hpp>

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
    Teuchos::Array<double> element_centroids( 
	d_mesh->dim() * d_mesh->localNumElements() );

    // Add the centroids from each block.
    int num_blocks = d_mesh->getNumBlocks();
    d_vertex_g2l.resize( num_blocks );
    Teuchos::RCP<Mesh> block;
    int block_num_elements = 0;
    int block_verts_per_elem = 0;
    Teuchos::ArrayRCP<GlobalOrdinal> block_vertex_ids;
    Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator vertex_id_it;
    int vlid = 0;
    Teuchos::ArrayRCP<GlobalOrdinal> block_vertex_coords;
    Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator vertex_coord_it;
    Teuchos::ArrayRCP<const GlobalOrdinal> block_connectivity;
    Teuchos::ArrayRCP<const GlobalOrdinal>::const_iterator connectivity_it;
    Teuchos::RCP<shards::CellTopology> block_topology;
    for ( int b = 0; b < num_blocks; ++i )
    {
	// Get the current block.
	block = mesh->getBlock( i );

	// Map the vertex ids in this block for faster coordinate access.
	block_vertex_ids = MeshTools<Mesh>::verticesNonConstView( *block );
	vlid = 0;
	for ( vertex_id_it = block_vertex_ids.begin();
	      vertex_id_it != block_vertex_ids.end();
	      ++vertex_id_it, ++vlid )
	{
	    d_vertex_g2l[b][*vertex_id_it] = vlid;
	}

	// Make a discretization cell for this block.
	block_topology = CellTopologyFactory::create(
	    MT::elementTopology(*block), MT::verticesPerElement(*block) );

	// Get the centroid of each element in the block.
	block_connectivity = MeshTools<Mesh>::connectivityView( *block );
	block_num_elements = MeshTools<Mesh>::numElements( *block );
	block_verts_per_elem = MT::verticesPerElement( *block );
	for ( int e = 0; e < block_num_elements; ++e )
	{
	    
	}
    }

    // Build a cloud search.
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
bool ElementTree<Mesh>::findPoint( const Teuchos::ArrayView<double>& coords,
				   GlobalOrdinal& element,
				   double tolerance )
{

}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ELEMENTTREE_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ElementTree_impl.hpp
//---------------------------------------------------------------------------//

