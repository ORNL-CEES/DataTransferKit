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
 * \file DTK_TopologyTools.cpp
 * \author Stuart R. Slattery
 * \brief TopologyTools definition.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "DTK_CellTopologyFactory.hpp"
#include "DTK_TopologyTools.hpp"
#include "DTK_GeometryTraits.hpp"
#include "DTK_DBC.hpp"
#include "DataTransferKit_config.hpp"

#include <MBGeomUtil.hpp>
#include <MBCartVect.hpp>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Box-element overlap query.
 *
 * \param box The box to check overlap with.
 *
 * \param element The element to check overlap with.
 *
 * \param moab The Moab interface containing the element.
 *
 * \return Return true if the box and element overlap, false if not.
 */
bool TopologyTools::boxElementOverlap( 
    const BoundingBox& box,
    const moab::EntityHandle element,
    const Teuchos::RCP<moab::Interface>& moab )
{
    // Get the element topology.
    moab::EntityType element_topology = moab->type_from_handle( element );

    // Get the element vertices.
    DTK_REMEMBER( moab::ErrorCode error );
    std::vector<moab::EntityHandle> element_vertices;
#if HAVE_DTK_DBC
    error = moab->get_adjacencies( &element,
				   1,
				   0,
				   false,
				   element_vertices );
#else
    moab->get_adjacencies( &element,
			   1,
			   0,
			   false,
			   element_vertices );
#endif
    DTK_CHECK( error == moab::MB_SUCCESS );

    // Extract the vertex coordinates.
    int num_element_vertices = element_vertices.size();
    Teuchos::Array<double> element_vertex_coords( 3 * num_element_vertices );
#if HAVE_DTK_DBC
    error = moab->get_coords( &element_vertices[0], 
			      num_element_vertices, 
			      &element_vertex_coords[0] );
#else
    moab->get_coords( &element_vertices[0], 
		      num_element_vertices, 
		      &element_vertex_coords[0] );
#endif
    DTK_CHECK( error == moab::MB_SUCCESS );

    // Extract box values.
    Teuchos::Tuple<double,6> bounds = box.getBounds();
    double half_dims[3] = { (bounds[3]-bounds[0])/2,
			    (bounds[4]-bounds[1])/2,
			    (bounds[5]-bounds[2])/2 };
    double center[3] = { half_dims[0] + bounds[0],
			 half_dims[1] + bounds[1],
			 half_dims[2] + bounds[2] };

    // Check for overlap.
    return moab::GeomUtil::box_elem_overlap( 
	(moab::CartVect*) &element_vertex_coords[0],
	element_topology,
	(moab::CartVect) center,
	(moab::CartVect) half_dims );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_TopologyTools.cpp
//---------------------------------------------------------------------------//
