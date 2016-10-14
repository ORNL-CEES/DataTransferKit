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
 * \file   DTK_IntrepidCell.cpp
 * \author Stuart Slattery
 * \brief  Manager for intrepid cell-level operations.
 */
//---------------------------------------------------------------------------//

#include "DTK_DBC.hpp"

#include "DTK_IntrepidCell.hpp"

#include <Teuchos_as.hpp>

#include <Intrepid_CellTools.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_Types.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
IntrepidCell::IntrepidCell( const shards::CellTopology &cell_topology,
                            const unsigned degree )
    : d_topology( cell_topology )
    , d_cub_points( 0, 0 )
    , d_cub_weights( 0 )
    , d_jacobian( 0, 0, 0, 0 )
    , d_jacobian_det( 0, 0 )
    , d_weighted_measures( 0, 0 )
    , d_physical_ip_coordinates( 0, 0, 0 )
{
    Intrepid::DefaultCubatureFactory<Scalar, MDArray> cub_factory;
    d_cubature = cub_factory.create( d_topology, degree );

    unsigned num_cub_points = d_cubature->getNumPoints();
    unsigned cub_dim = cell_topology.getDimension();

    d_cub_points.resize( num_cub_points, cub_dim );
    d_cub_weights.resize( num_cub_points );
    d_cubature->getCubature( d_cub_points, d_cub_weights );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
IntrepidCell::~IntrepidCell() { /* ... */}

//---------------------------------------------------------------------------//
/*!
 * \brief Given physical coordinates for the cell nodes (Cell,Node,Dim),
 * assign them to the cell without allocating internal data.
 */
void IntrepidCell::setCellNodeCoordinates( const MDArray &cell_node_coords )
{
    d_cell_node_coords = cell_node_coords;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given physical coordinates for the cell nodes (Cell,Node,Dim),
 * allocate the state of the cell object.
 */
void IntrepidCell::allocateCellState( const MDArray &cell_node_coords )
{
    // Store the cell node coords as the current state.
    setCellNodeCoordinates( cell_node_coords );

    // Get required dimensions.
    int num_cells = d_cell_node_coords.dimension( 0 );
    int num_ip = d_cub_points.dimension( 0 );
    int space_dim = d_cub_points.dimension( 1 );

    // Resize arrays.
    d_jacobian.resize( num_cells, num_ip, space_dim, space_dim );
    d_jacobian_det.resize( num_cells, num_ip );
    d_weighted_measures.resize( num_cells, num_ip );
    d_physical_ip_coordinates.resize( num_cells, num_ip, space_dim );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the state of the cell object for the current cell node
 * coordinates.
 */
void IntrepidCell::updateCellState()
{
    // Compute the Jacobian.
    Intrepid::CellTools<Scalar>::setJacobian( d_jacobian, d_cub_points,
                                              d_cell_node_coords, d_topology );

    // Compute the determinant of the Jacobian.
    Intrepid::CellTools<Scalar>::setJacobianDet( d_jacobian_det, d_jacobian );

    // Compute the cell measures.
    Intrepid::FunctionSpaceTools::computeCellMeasure<Scalar>(
        d_weighted_measures, d_jacobian_det, d_cub_weights );

    // Compute physical frame integration point coordinates.
    Intrepid::CellTools<Scalar>::mapToPhysicalFrame(
        d_physical_ip_coordinates, d_cub_points, d_cell_node_coords,
        d_topology );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Free function for updating the cell state for a new set of
 * physical cells in a single call.
 */
void IntrepidCell::updateState( IntrepidCell &intrepid_cell,
                                const MDArray &cell_node_coords )
{
    intrepid_cell.allocateCellState( cell_node_coords );
    intrepid_cell.updateCellState();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a set of coordinates in the physical frame of the cell, map
 * them to the reference frame of the cell.
 */
void IntrepidCell::mapToCellReferenceFrame( const MDArray &physical_coords,
                                            MDArray &reference_coords )
{
    DTK_REQUIRE( 2 == physical_coords.rank() );
    DTK_REQUIRE( 2 == reference_coords.rank() );
    Intrepid::CellTools<Scalar>::mapToReferenceFrame(
        reference_coords, physical_coords, d_cell_node_coords, d_topology, 0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a set of coordinates in the reference frame of the cell, map
 * them to the physical frame.
 */
void IntrepidCell::mapToCellPhysicalFrame( const MDArray &parametric_coords,
                                           MDArray &physical_coords )
{
    DTK_REQUIRE( 2 == parametric_coords.rank() );
    DTK_REQUIRE( 3 == physical_coords.rank() );
    DTK_REQUIRE( parametric_coords.dimension( 1 ) ==
                 Teuchos::as<int>( d_topology.getDimension() ) );
    DTK_REQUIRE( physical_coords.dimension( 0 ) ==
                 d_cell_node_coords.dimension( 0 ) );
    DTK_REQUIRE( physical_coords.dimension( 1 ) ==
                 parametric_coords.dimension( 0 ) );
    DTK_REQUIRE( physical_coords.dimension( 2 ) ==
                 Teuchos::as<int>( d_topology.getDimension() ) );

    Intrepid::CellTools<Scalar>::mapToPhysicalFrame(
        physical_coords, parametric_coords, d_cell_node_coords, d_topology );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point given in parametric coordinates is inside of the
 * reference cell.
 */
bool IntrepidCell::pointInReferenceCell( const MDArray &reference_point,
                                         const double tolerance )
{
    return Intrepid::CellTools<Scalar>::checkPointsetInclusion(
        reference_point, d_topology, tolerance );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point given in physical coordinates is inside of the
 * phyiscal cell.
 */
bool IntrepidCell::pointInPhysicalCell( const MDArray &point,
                                        const double tolerance )
{
    MDArray reference_point( point );
    reference_point.initialize();
    mapToCellReferenceFrame( point, reference_point );
    return pointInReferenceCell( reference_point, tolerance );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of cells in the current state.
 */
int IntrepidCell::getNumCells() const
{
    return d_weighted_measures.dimension( 0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of cell points.
 */
int IntrepidCell::getNumIntegrationPoints() const
{
    return d_cub_points.dimension( 0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the spatial dimension.
 */
int IntrepidCell::getSpatialDimension() const
{
    return d_cub_points.dimension( 1 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the cell measures (Cell). cell_measures must all ready be
 * allocated.
 */
void IntrepidCell::getCellMeasures( MDArray &cell_measures ) const
{
    DTK_REQUIRE( 1 == cell_measures.rank() );
    DTK_REQUIRE( cell_measures.dimension( 0 ) ==
                 d_weighted_measures.dimension( 0 ) );

    MDArray dofs( d_cell_node_coords.dimension( 0 ),
                  d_cub_weights.dimension( 0 ) );
    dofs.initialize( 1.0 );
    integrate( dofs, cell_measures );
}

//---------------------------------------------------------------------------//
// Get the physical cell point coordinates in each cell
// (Cell,IP,Dim).
void IntrepidCell::getPhysicalIntegrationCoordinates(
    MDArray &physical_ip_coordinates ) const
{
    physical_ip_coordinates = d_physical_ip_coordinates;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given DOFs at the quadrature points {(Cell,Node) for scalar fields,
 * (Cell,Node,VecDim) for vector fields, and (Cell,Node,TensDim1,TensDim2) for
 * tensor fields.} perform the numerical integration in each cell by
 * contracting them with the weighted measures.
 */
void IntrepidCell::integrate( const MDArray &dofs, MDArray &integrals ) const
{
    Intrepid::FunctionSpaceTools::integrate<Scalar>(
        integrals, dofs, d_weighted_measures, Intrepid::COMP_BLAS );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_IntrepidCell.cpp
//---------------------------------------------------------------------------//
