//---------------------------------------------------------------------------//
/*
  Copyright (c) 2013, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \file   DTK_IntrepidSideCell_impl.hpp
 * \author Stuart Slattery
 * \brief  Manager for intrepid cell-level operations on cell sides.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTREPIDSIDECELL_IMPL_HPP
#define DTK_INTREPIDSIDECELL_IMPL_HPP

#include "DTK_DBC.hpp"

#include <Teuchos_as.hpp>

#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename MDArray>
IntrepidSideCell<MDArray>::IntrepidSideCell( 
    const shards::CellTopology& side_topology,
    const unsigned side_id,
    const shards::CellTopology& parent_topology,
    const unsigned degree )
    : Base( side_topology, degree )
    , d_side_id( side_id )
    , d_parent_topology( parent_topology )
{
    // Map the side cubature points to the cell frame.
    MDArray mapped_cub_points( this->d_cubature->getNumPoints(), 
			       parent_topology.getDimension() );
    Intrepid::CellTools<Scalar>::mapToReferenceSubcell(
	mapped_cub_points, this->d_cub_points, side_topology.getDimension(),
	d_side_id, d_parent_topology );
    this->d_cub_points = mapped_cub_points;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename MDArray>
IntrepidSideCell<MDArray>::~IntrepidSideCell()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Update the cell state of the object for the current cell node
 * coordinates.
 */
template<typename MDArray>
void IntrepidSideCell<MDArray>::updateCellState()
{
    unsigned space_dim = d_parent_topology.getDimension();

    // Compute the Jacobian.
    Intrepid::CellTools<Scalar>::setJacobian( 
	this->d_jacobian, this->d_cub_points, 
	this->d_cell_node_coords, d_parent_topology );

    // Compute the cell side measures.
    switch ( space_dim )
    {

	case 3:
	{
	    // Face case.
	    Intrepid::FunctionSpaceTools::computeFaceMeasure<Scalar>( 
		this->d_weighted_measures, this->d_jacobian, 
		this->d_cub_weights, d_side_id, d_parent_topology );
	}
	break;

	case 2:
	{
	    // Edge case.
	    Intrepid::FunctionSpaceTools::computeEdgeMeasure<Scalar>( 
		this->d_weighted_measures, this->d_jacobian,
		this->d_cub_weights, d_side_id, d_parent_topology );
	}
	break;

	default:
	    DTK_INSIST( 3 == space_dim || 2 == space_dim,
			  "Subcell cell not supported for dimension" );
	    break;
    }

    // Compute physical frame integration point coordinates.
    Intrepid::CellTools<Scalar>::mapToPhysicalFrame( 
	this->d_physical_ip_coordinates, this->d_cub_points, 
	this->d_cell_node_coords, d_parent_topology );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a set of coordinates in the reference frame of the cell, map
 * them to the physical frame.
 */
template<typename MDArray>
void IntrepidSideCell<MDArray>::mapToCellPhysicalFrame( 
    const MDArray& parametric_coords, MDArray& physical_coords )
{
    DTK_REQUIRE( 2 == parametric_coords.rank() );
    DTK_REQUIRE( 3 == physical_coords.rank() );
    DTK_REQUIRE( parametric_coords.dimension(1) ==
		   Teuchos::as<int>(this->d_topology.getDimension()) );
    DTK_REQUIRE( physical_coords.dimension(0) ==
		   this->d_cell_node_coords.dimension(0) );
    DTK_REQUIRE( physical_coords.dimension(1) ==
		   parametric_coords.dimension(0) );
    DTK_REQUIRE( physical_coords.dimension(2) ==
		   Teuchos::as<int>(d_parent_topology.getDimension()) );

    MDArray mapped_coords( parametric_coords.dimension(0),
			   d_parent_topology.getDimension() );

    Intrepid::CellTools<Scalar>::mapToReferenceSubcell(
	mapped_coords, parametric_coords, this->d_topology.getDimension(),
	d_side_id, d_parent_topology );

    Intrepid::CellTools<Scalar>::mapToPhysicalFrame( 
	physical_coords, mapped_coords, 
	this->d_cell_node_coords, d_parent_topology );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the physical normals of the side.
 */
template<typename MDArray>
void IntrepidSideCell<MDArray>::getPhysicalSideNormals( MDArray& side_normals )
{
    Intrepid::CellTools<Scalar>::getPhysicalSideNormals(
	side_normals, this->d_jacobian, d_side_id, d_parent_topology );
}

//---------------------------------------------------------------------------//

} // End namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTREPIDSIDECELL_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_IntrepidSideCell.hpp
//---------------------------------------------------------------------------//

