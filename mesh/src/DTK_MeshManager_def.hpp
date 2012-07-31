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

#include "DTK_MeshTypes.hpp"
#include "DTK_MeshTools.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Mesh>
MeshManager<Mesh>::MeshManager( const Teuchos::ArrayRCP<Mesh>& mesh_blocks,
				const RCP_Comm& comm,
				const std::size_t dim )
    : d_mesh_blocks( mesh_blocks )
    , d_comm( comm )
    , d_dim( dim )
    , d_active_nodes( d_mesh_blocks.size() )
    , d_active_elements( d_mesh_blocks.size() )
{
    validate();
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
 * \brief Get the local number of elements in the mesh.
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
	local_num_elements += MeshTools<Mesh>::numElements( *block_iterator );
    }
    return local_num_elements;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global number of elements in the mesh.
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
 * \brief Compute the global bounding box around the entire mesh.
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
	if ( MeshTools<Mesh>::numNodes( *block_iterator ) > 0 )
	{
	    block_box =
		MeshTools<Mesh>::globalBoundingBox( *block_iterator, d_comm );

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
    BlockIterator block_iterator;
    for ( block_iterator = d_mesh_blocks.begin();
	  block_iterator != d_mesh_blocks.end();
	  ++block_iterator )
    {
	testInvariant( d_dim == MT::nodeDim( *block_iterator ),
		       "Mesh dimension != node dimension" );

	if ( d_dim == 0 )
	{
	    testInvariant( MT::elementTopology( *block_iterator ) == DTK_VERTEX,
			   "Element topology does not match mesh dimension" );
	}
	else if ( d_dim == 1 )
	{
	    testInvariant( MT::elementTopology( *block_iterator ) == 
			   DTK_LINE_SEGMENT,
			   "Element topology does not match mesh dimension" );
	}
	else if ( d_dim == 2 )
	{
	    testInvariant( MT::elementTopology( *block_iterator ) == 
			   DTK_TRIANGLE ||
			   MT::elementTopology( *block_iterator ) == 
			   DTK_QUADRILATERAL,
			   "Element topology does not match mesh dimension" );
	}
	else if ( d_dim == 3 )
	{
	    testInvariant( MT::elementTopology( *block_iterator ) == 
			   DTK_TETRAHEDRON ||
			   MT::elementTopology( *block_iterator ) == 
			   DTK_HEXAHEDRON ||
			   MT::elementTopology( *block_iterator ) == 
			   DTK_PYRAMID ||
			   MT::elementTopology( *block_iterator ) == 
			   DTK_WEDGE,
			   "Element topology does not match mesh dimension" );
	}
	else
	{
	    throw MeshException( "Mesh dimension > 3 not supported" );
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESHMANAGER_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshManager_def.hpp
//---------------------------------------------------------------------------//


