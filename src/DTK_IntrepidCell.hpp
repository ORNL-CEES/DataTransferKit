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
 * \file   DTK_IntrepidCell.hpp
 * \author Stuart Slattery
 * \brief  Manager for Intrepid cell-level operations.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTREPIDCELL_HPP
#define DTK_INTREPIDCELL_HPP

#include "DTK_Tolerances.hpp"

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_Cubature.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class IntrepidCell
 * \brief Manager for Intrepid cell-level operations.
 */
//---------------------------------------------------------------------------//
template<typename MDArray>
class IntrepidCell
{
  public:

    //@{
    //! Typedefs.
    typedef typename MDArray::scalar_type Scalar;
    //@}

  public:

    // Constructor.
    IntrepidCell( const shards::CellTopology& cell_topology,
		  const unsigned degree );

    // Destructor.
    virtual ~IntrepidCell();

    // Given physical coordinates for the cell nodes (Cell,Node,Dim), assign
    // them to the cell without allocating internal data.
    void setCellNodeCoordinates( const MDArray& cell_node_coords );

    // Given physical coordinates for the cell nodes (Cell,Node,Dim),
    // allocate the state of the cell object.
    void allocateCellState( const MDArray& cell_node_coords );

    // Update the cell state of the object for the current cell node
    // coordinates.
    virtual void updateCellState();

    // Free function for updating the cell state for a new set of
    // physical cells in a single call.
    static void updateState( IntrepidCell& Intrepid_cell,
			     const MDArray& cell_node_coords );

    // Given a set of coordinates in the reference frame of the cell
    // (Node,Dim), map them to the physical frame (Cell,Node,Dim).
    virtual void mapToCellPhysicalFrame( const MDArray& parametric_coords,
					 MDArray& physical_coords );

    // Determine if a point given in natural coordinates is inside of the
    // reference cell.
    bool pointInReferenceCell( const MDArray& reference_point,
			       const double tolerance );

    // Determine if a point in physical coordinates is inside of the physical
    // cell.
    bool pointInPhysicalCell( const MDArray& point,
			      const double tolerance );

    // Get the number of cells in the current state.
    int getNumCells() const;

    // Get the number of cell points.
    int getNumIntegrationPoints() const;

    // Get the spatial dimension.
    int getSpatialDimension() const;

    // Get the cell measures (Cell). cell_measures must already be allocated.
    void getCellMeasures( MDArray& cell_measures ) const;

    // Get the physical integration point coordinates in each cell
    // (Cell,IP,Dim).
    void getPhysicalIntegrationCoordinates(
	MDArray& physical_ip_coordinates ) const;

    // Given DOFs at the quadrature points {(Cell,Node) for scalar fields,
    // (Cell,Node,VecDim) for vector fields, and (Cell,Node,TensDim1,TensDim2)
    // for tensor fields.} perform the numerical integration in each cell by
    // contracting them with the weighted measures.
    void integrate( const MDArray& dofs, MDArray& integrals ) const;

  protected:

    // Cell topology.
    shards::CellTopology d_topology;

    // Cell integration rule.
    Teuchos::RCP<Intrepid::Cubature<Scalar,MDArray> > d_cubature;

    // Cubature points (IP,Dim).
    MDArray d_cub_points;

    // Cubature weights (IP).
    MDArray d_cub_weights;

    // Jacobian (Cell,IP,Dim,Dim).
    MDArray d_jacobian;

    // Jacobian inverse (Cell,IP,Dim,Dim).
    MDArray d_jacobian_inv;

    // Jacobian determinant (Cell,IP).
    MDArray d_jacobian_det;

    // Weighted cell measures (Cell,IP).
    MDArray d_weighted_measures;

    // Current physical cell node coordinates.
    MDArray d_cell_node_coords;

    // Current physical integration point coordinates.
    MDArray d_physical_ip_coordinates;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_IntrepidCell_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_INTREPIDCELL_HPP

//---------------------------------------------------------------------------//
// end DTK_IntrepidCell.hpp
//---------------------------------------------------------------------------//

