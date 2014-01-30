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
 * \file DTK_TopologyTools_def.hpp
 * \author Stuart R. Slattery
 * \brief TopologyTools template definitions.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TOPOLOGYTOOLS_DEF_HPP
#define DTK_TOPOLOGYTOOLS_DEF_HPP

#include "DTK_GeometryTraits.hpp"
#include "DTK_DBC.hpp"
#include "DataTransferKit_config.hpp"

#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Element-in-geometry query.
 *
 * \param geometry The geometry.
 *
 * \param element The element. 
 *
 * \param moab The Moab interface containing the element.
 *
 * \param tolerance Tolerance used for element vertex-in-geometry
 * checks.
 *
 * \param all_vertices_for_inclusion Flag for element-in-geometry
 * inclusion. If set to true, all of an element's vertices are required to
 * reside within a geometry within the geometric tolerance in order to be
 * considered a member of that geometry's conformal mesh. If set to false,
 * only one of an element's vertices must be contained within the geometric
 * tolerance of the geometry in order to be considered a member of that
 * geometry's conformal mesh.
 *
 * \return Return true if any of the element's vertices are in the
 * geometry. This is based on the conformal mesh/geometry assumption.
 */
template<class Geometry>
bool TopologyTools::elementInGeometry( 
    const Geometry& geometry,
    const moab::EntityHandle element,
    const Teuchos::RCP<moab::Interface>& moab,
    const double tolerance,
    bool all_vertices_for_inclusion )
{
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

    // Check the vertex coordinates for inclusion in the geometry. 
    Teuchos::Array<double> vertex_coords(3);
    int verts_in_geometry = 0;
    for ( int i = 0; i < num_element_vertices; ++i )
    {
	vertex_coords[0] = element_vertex_coords[3*i];
	vertex_coords[1] = element_vertex_coords[3*i + 1];
	vertex_coords[2] = element_vertex_coords[3*i + 2];

	if ( GeometryTraits<Geometry>::pointInGeometry( geometry, 
							vertex_coords,
							tolerance ) )
	{
	    // All vertices required for inclusion case.
	    if ( all_vertices_for_inclusion )
	    {
		++verts_in_geometry;
	    }
	    // Only one vertex required for inclusion case.
	    else
	    {
		return true;
	    }
	}
    }

    if ( verts_in_geometry == num_element_vertices )
    { 
	return true;
    }
    else
    {
	return false;
    }
}

//---------------------------------------------------------------------------//

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_TOPOLOGYTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_TopologyTools_def.hpp
//---------------------------------------------------------------------------//
