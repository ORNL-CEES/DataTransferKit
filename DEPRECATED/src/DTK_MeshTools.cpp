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
 * \file DTK_MeshTools.cpp
 * \author Stuart R. Slattery
 * \brief MeshTools definition.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <iterator>

#include "DTK_DBC.hpp"
#include "DTK_MeshTools.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box for a mesh block.
 * 
 * \param mesh The mesh block for which to get the local bounding box.
 *
 * \return The local bounding box around the mesh block.
 */
BoundingBox MeshTools::localBoundingBox( const Teuchos::RCP<MeshBlock>& mesh )
{
    double x_min = -Teuchos::ScalarTraits<double>::rmax();
    double y_min = -Teuchos::ScalarTraits<double>::rmax();
    double z_min = -Teuchos::ScalarTraits<double>::rmax();

    double x_max = Teuchos::ScalarTraits<double>::rmax();
    double y_max = Teuchos::ScalarTraits<double>::rmax();
    double z_max = Teuchos::ScalarTraits<double>::rmax();

    MeshId num_vertices = mesh->numVertices();
    int vertex_dim = mesh->dimension();

    if ( vertex_dim > 0 )
    {
	x_min = *std::min_element( mesh->vertexCoordinates().begin(),
				   mesh->vertexCoordinates().begin() + num_vertices );
	x_max = *std::max_element( 
	    mesh->vertexCoordinates().begin(),
	    mesh->vertexCoordinates().begin() + num_vertices );
    }
    if ( vertex_dim > 1 )
    {
	y_min = *std::min_element( 
	    mesh->vertexCoordinates().begin() + num_vertices,
	    mesh->vertexCoordinates().begin() + 2*num_vertices );
	y_max = *std::max_element( 
	    mesh->vertexCoordinates().begin() + num_vertices,
	    mesh->vertexCoordinates().begin() + 2*num_vertices );
    }
    if ( vertex_dim > 2 )
    {
	z_min = *std::min_element( 
	    mesh->vertexCoordinates().begin() + 2*num_vertices,
	    mesh->vertexCoordinates().begin() + 3*num_vertices );
	z_max = *std::max_element( 
	    mesh->vertexCoordinates().begin() + 2*num_vertices,
	    mesh->vertexCoordinates().begin() + 3*num_vertices );
    }

    return BoundingBox( x_min, y_min, z_min, x_max, y_max, z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global bounding box for a mesh block over the given
 * communicator.
 *
 * \param mesh The mesh block over which to get the global bounding box.
 *
 * \param comm The communicator over which the mesh block is defined.
 *
 * \return The global bounding box over the mesh block.
 */
BoundingBox MeshTools::globalBoundingBox( const Teuchos::RCP<MeshBlock>& mesh, 
					  const RCP_Comm& comm )
{
    BoundingBox local_box = localBoundingBox( mesh );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    double min_bounds[3];
    double max_bounds[3];

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    3,
				    &local_bounds[0],
				    min_bounds );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    3,
				    &local_bounds[3],
				    max_bounds );

    return BoundingBox( min_bounds[0], min_bounds[1], min_bounds[2],
			max_bounds[0], max_bounds[1], max_bounds[2] );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MeshTools.cpp
//---------------------------------------------------------------------------//

