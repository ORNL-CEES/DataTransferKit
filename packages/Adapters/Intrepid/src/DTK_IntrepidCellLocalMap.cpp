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
 * \brief DTK_IntrepidCellLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Helper functions for implementing the local map interface with
 * intrepid.
 */
//---------------------------------------------------------------------------//

#include "DTK_IntrepidCellLocalMap.hpp"
#include "DTK_DBC.hpp"
#include "DTK_IntrepidCell.hpp"
#include "DTK_ProjectionPrimitives.hpp"

#include <Intrepid_FieldContainer.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
double IntrepidCellLocalMap::measure(
    const shards::CellTopology &entity_topo,
    const Intrepid::FieldContainer<double> &entity_coords )
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell entity_cell( entity_topo, 1 );

    // Update thet state of the cell.
    IntrepidCell::updateState( entity_cell, entity_coords );

    // Compute the measure of the cell.
    Intrepid::FieldContainer<double> measure( 1 );
    entity_cell.getCellMeasures( measure );
    return measure( 0 );
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void IntrepidCellLocalMap::centroid(
    const shards::CellTopology &entity_topo,
    const Intrepid::FieldContainer<double> &entity_coords,
    const Teuchos::ArrayView<double> &centroid )
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell entity_cell( entity_topo, 1 );
    entity_cell.setCellNodeCoordinates( entity_coords );

    // Get the reference center of the cell.
    Intrepid::FieldContainer<double> ref_center( 1,
                                                 entity_coords.dimension( 2 ) );
    ProjectionPrimitives::referenceCellCenter( entity_topo, ref_center );

    // Map the cell center to the physical frame.
    Intrepid::FieldContainer<double> phys_center(
        1, 1, entity_coords.dimension( 2 ) );
    entity_cell.mapToCellPhysicalFrame( ref_center, phys_center );

    // Extract the centroid coordinates.
    centroid.assign( phys_center.getData()() );
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
bool IntrepidCellLocalMap::mapToReferenceFrame(
    const shards::CellTopology &entity_topo,
    const Intrepid::FieldContainer<double> &entity_coords,
    const Teuchos::ArrayView<const double> &point,
    const Teuchos::ArrayView<double> &reference_point )
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell entity_cell( entity_topo, 1 );
    entity_cell.setCellNodeCoordinates( entity_coords );

    // Map the point to the reference frame of the cell.
    Teuchos::Array<int> array_dims( 2 );
    array_dims[0] = 1;
    array_dims[1] = entity_coords.dimension( 2 );
    Intrepid::FieldContainer<double> point_container(
        array_dims, const_cast<double *>( point.getRawPtr() ) );
    Intrepid::FieldContainer<double> ref_point_container(
        array_dims, reference_point.getRawPtr() );
    entity_cell.mapToCellReferenceFrame( point_container, ref_point_container );

    // Return true to indicate successful mapping. Catching Intrepid errors
    // and returning false is a possibility here.
    return true;
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool IntrepidCellLocalMap::checkPointInclusion(
    const shards::CellTopology &entity_topo,
    const Teuchos::ArrayView<const double> &reference_point,
    const double tolerance )
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell entity_cell( entity_topo, 1 );

    // Check point inclusion.
    Teuchos::Array<int> array_dims( 2 );
    array_dims[0] = 1;
    array_dims[1] = reference_point.size();
    Intrepid::FieldContainer<double> ref_point_container(
        array_dims, const_cast<double *>( reference_point.getRawPtr() ) );
    return entity_cell.pointInReferenceCell( ref_point_container, tolerance );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void IntrepidCellLocalMap::mapToPhysicalFrame(
    const shards::CellTopology &entity_topo,
    const Intrepid::FieldContainer<double> &entity_coords,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &point )
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell entity_cell( entity_topo, 1 );
    entity_cell.setCellNodeCoordinates( entity_coords );

    // Map the reference point to the physical frame of the cell.
    Teuchos::Array<int> ref_array_dims( 2 );
    ref_array_dims[0] = 1;
    ref_array_dims[1] = entity_coords.dimension( 2 );
    Intrepid::FieldContainer<double> ref_point_container(
        ref_array_dims, const_cast<double *>( reference_point.getRawPtr() ) );
    Teuchos::Array<int> phys_array_dims( 3 );
    phys_array_dims[0] = 1;
    phys_array_dims[1] = 1;
    phys_array_dims[2] = entity_coords.dimension( 2 );
    Intrepid::FieldContainer<double> point_container( phys_array_dims,
                                                      point.getRawPtr() );
    entity_cell.mapToCellPhysicalFrame( ref_point_container, point_container );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_IntrepidCellLocalMap.cpp
//---------------------------------------------------------------------------//
