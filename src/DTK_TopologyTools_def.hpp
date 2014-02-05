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

#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get the center of the reference cell of the given topology.
 */
template<typename MDArray>
void TopologyTools::referenceCellCenter( 
    const shards::CellTopology& cell_topo, MDArray& cell_center )
{
    typedef typename MDArray::scalar_type Scalar;

    DTK_REQUIRE( 2 == cell_center.rank() );
    DTK_REQUIRE( Teuchos::as<unsigned>(cell_center.dimension(1)) == 
		   cell_topo.getDimension() );

    int num_cells = cell_center.dimension(0);    

    switch( cell_topo.getKey() ){
	case shards::Line<2>::key:
	case shards::Line<3>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = Teuchos::ScalarTraits<Scalar>::zero();
	    }
	    break;
      
	case shards::Triangle<3>::key:
	case shards::Triangle<4>::key:
	case shards::Triangle<6>::key:    
	    for ( int n = 0; n < num_cells; ++n )
	    {

		cell_center(n,0) = 1.0/3.0;
		cell_center(n,1) = 1.0/3.0;  
	    }
	    break;
      
	case shards::Quadrilateral<4>::key:
	case shards::Quadrilateral<8>::key:
	case shards::Quadrilateral<9>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = Teuchos::ScalarTraits<Scalar>::zero();      
		cell_center(n,1) = Teuchos::ScalarTraits<Scalar>::zero();    
	    }
	    break;
      
	case shards::Tetrahedron<4>::key:
	case shards::Tetrahedron<10>::key:
	case shards::Tetrahedron<11>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = 1.0/6.0;    
		cell_center(n,1) = 1.0/6.0;    
		cell_center(n,2) = 1.0/6.0; 
	    }
	    break;
      
	case shards::Hexahedron<8>::key:
	case shards::Hexahedron<20>::key:
	case shards::Hexahedron<27>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = Teuchos::ScalarTraits<Scalar>::zero();
		cell_center(n,1) = Teuchos::ScalarTraits<Scalar>::zero();
		cell_center(n,2) = Teuchos::ScalarTraits<Scalar>::zero();
	    }
	    break;
      
	case shards::Wedge<6>::key:
	case shards::Wedge<15>::key:
	case shards::Wedge<18>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = 1.0/3.0;
		cell_center(n,1) = 1.0/3.0;
		cell_center(n,2) = Teuchos::ScalarTraits<Scalar>::zero();
	    }
	    break;

	default:
	    DTK_INSIST( false );
	    break;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-in-geometry query.
 *
 * \param geometry The geometry.
 *
 * \param element_node_coords The coordinates of the element vertices.
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
template<class Geometry, typename MDArray>
bool TopologyTools::elementInGeometry( 
    const Geometry& geometry,
    const MDArray& element_node_coords,
    const double tolerance,
    bool all_vertices_for_inclusion )
{
    int space_dim = element_node_coords.dimension(0);
    int num_element_vertices = element_node_coords.dimension(1);
    Teuchos::Array<double> vertex_coords( space_dim );
    int verts_in_geometry = 0;
    for ( int n = 0; n < num_element_vertices; ++n )
    {
	for ( int d = 0; d < space_dim; ++d )
	{
	    vertex_coords[d] = element_node_coords(0,n,d);
	}

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
    
    return false;
}

//---------------------------------------------------------------------------//

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_TOPOLOGYTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_TopologyTools_def.hpp
//---------------------------------------------------------------------------//
