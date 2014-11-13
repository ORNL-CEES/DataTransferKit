//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
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
 * \brief DTK_EntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshEntityLocalMap.hpp"
#include "DTK_STKMeshHelpers.hpp"
#include "DTK_IntrepidSideCell.hpp"
#include "DTK_DBC.hpp"

#include <Intrepid_FieldContainer.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
STKMeshEntityLocalMap::STKMeshEntityLocalMap(
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data )
    : d_bulk_data( bulk_data )
{ 
    // Get the mesh parts.
    const stk::mesh::PartVector& all_parts =
	d_bulk_data->mesh_meta_data().get_parts();

    // Construct the Intrepid cells for the parts that are not empty.
    for ( auto part_it = all_parts.begin(); 
	  part_it != all_parts.end(); 
	  ++part_it )
    {
	// Make a selector from the part.
	stk::mesh::Selector part_select( **part_it );

	// If the part is not empty, make a cell for the part.
	if ( !part_select.is_empty(stk::topology::ELEM_RANK) )
	{
	    stk::topology stk_part_topo = (*part_it)->topology();
	    if ( stk::topology::INVALID_TOPOLOGY != stk_part_topo )
	    {
		shards::CellTopology part_topo = 
		    stk::mesh::get_cell_topology( stk_part_topo );
		if ( !d_topo_to_cell_map.count(part_topo.getKey()) )
		{
		    d_intrepid_cells.push_back( 
			Teuchos::rcp(new IntrepidCell(part_topo,1)) );
		    d_topo_to_cell_map.insert( 
			std::make_pair(part_topo.getKey(),
				       d_intrepid_cells.size()-1) );
		}
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Destructor.
STKMeshEntityLocalMap::~STKMeshEntityLocalMap()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D). 
double STKMeshEntityLocalMap::measure( const Entity& entity ) const
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell& entity_cell = getEntityIntrepidCell( entity );

    // Load the node coordinates of the entity into the cell.
    Intrepid::FieldContainer<double> entity_coords = 
	STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(1,STKMeshHelpers::extractEntity(entity)),
	    *d_bulk_data );
    IntrepidCell::updateState( entity_cell, entity_coords );
    
    // Compute the measure of the cell.
    Intrepid::FieldContainer<double> measure(1);
    entity_cell.getCellMeasures( measure );
    return measure(0);
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void STKMeshEntityLocalMap::centroid( 
    const Entity& entity, const Teuchos::ArrayView<double>& centroid ) const
{ 
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell& entity_cell = getEntityIntrepidCell( entity );

    // Load the node coordinates of the entity into the cell.
    Intrepid::FieldContainer<double> entity_coords = 
	STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(1,STKMeshHelpers::extractEntity(entity)),
	    *d_bulk_data );
    entity_cell.setCellNodeCoordinates( entity_coords );

    // Get the reference center of the cell.
    Intrepid::FieldContainer<double> ref_center(1,entity.physicalDimension());
    referenceCellCenter( entity, ref_center );

    // Map the cell center to the physical frame.
    Intrepid::FieldContainer<double> phys_center(1,1,entity.physicalDimension());
    entity_cell.mapToCellPhysicalFrame( ref_center, phys_center );
    
    // Extract the centroid coordinates.
    centroid.assign( phys_center.getData()() );
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance. 
bool STKMeshEntityLocalMap::isSafeToMapToReferenceFrame(
    const Entity& entity,
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::RCP<MappingStatus>& status ) const
{
    // Get the bounding box of the entity.
    Teuchos::Tuple<double,6> entity_box;
    entity.boundingBox( entity_box );

    // Get the test tolerance.
    double tolerance = 1.0e-6;
    if ( Teuchos::nonnull(this->b_parameters) )
    {
	if ( this->b_parameters->isParameter("Point Inclusion Tolerance") )
	{
	    tolerance = 
		this->b_parameters->get<double>("Point Inclusion Tolerance");
	}
    }

    // Check if the point is in the bounding box of the entity.
    int space_dim = entity.physicalDimension();
    bool in_x = true;
    if ( space_dim > 0 )
    {
	double x_tol = (entity_box[3] - entity_box[0])*tolerance;
	in_x = ( (point[0] >= (entity_box[0] - x_tol)) &&
		 (point[0] <= (entity_box[3] + x_tol)) );
    }
    bool in_y = true;
    if ( space_dim > 1 )
    {
	double y_tol = (entity_box[4] - entity_box[1])*tolerance;
	in_y = ( (point[1] >= (entity_box[1] - y_tol)) &&
		 (point[1] <= (entity_box[4] + y_tol)) );
    }
    bool in_z = true;
    if ( space_dim > 2 )
    {
	double z_tol = (entity_box[5] - entity_box[2])*tolerance;
	in_z = ( (point[2] >= (entity_box[2] - z_tol)) &&
		 (point[2] <= (entity_box[5] + z_tol)) );
    }
    return (in_x && in_y && in_z);
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
bool STKMeshEntityLocalMap::mapToReferenceFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point,
    const Teuchos::RCP<MappingStatus>& status ) const
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell& entity_cell = getEntityIntrepidCell( entity );

    // Load the node coordinates of the entity into the cell.
    Intrepid::FieldContainer<double> entity_coords = 
	STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(1,STKMeshHelpers::extractEntity(entity)),
	    *d_bulk_data );
    entity_cell.setCellNodeCoordinates( entity_coords );

    // Map the point to the reference frame of the cell.
    Teuchos::Array<int> array_dims(2);
    array_dims[0] = 1;
    array_dims[1] = entity.physicalDimension();
    Intrepid::FieldContainer<double> point_container( 
	array_dims, const_cast<double*>(point.getRawPtr()) );
    Intrepid::FieldContainer<double> ref_point_container( 
	array_dims, reference_point.getRawPtr() );
    entity_cell.mapToCellReferenceFrame( point_container, ref_point_container );

    // Return true to indicate successful mapping. Catching Intrepid errors
    // and returning false is a possibility here.
    return true;
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool STKMeshEntityLocalMap::checkPointInclusion( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    // Get the test tolerance.
    double tolerance = 1.0e-6;
    if ( Teuchos::nonnull(this->b_parameters) )
    {
	if ( this->b_parameters->isParameter("Point Inclusion Tolerance") )
	{
	    tolerance = 
		this->b_parameters->get<double>("Point Inclusion Tolerance");
	}
    }

    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell& entity_cell = getEntityIntrepidCell( entity );

    // Check point inclusion.
    Teuchos::Array<int> array_dims(2);
    array_dims[0] = 1;
    array_dims[1] = entity.physicalDimension();
    Intrepid::FieldContainer<double> ref_point_container( 
	array_dims, const_cast<double*>(reference_point.getRawPtr()) );
    return entity_cell.pointInReferenceCell( ref_point_container, tolerance );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void STKMeshEntityLocalMap::mapToPhysicalFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    // Get the Intrepid cell corresponding to the entity topology.
    IntrepidCell& entity_cell = getEntityIntrepidCell( entity );

    // Load the node coordinates of the entity into the cell.
    Intrepid::FieldContainer<double> entity_coords = 
	STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(1,STKMeshHelpers::extractEntity(entity)),
	    *d_bulk_data );
    entity_cell.setCellNodeCoordinates( entity_coords );

    // Map the reference point to the physical frame of the cell.
    Teuchos::Array<int> ref_array_dims(2);
    ref_array_dims[0] = 1;
    ref_array_dims[1] = entity.physicalDimension();
    Intrepid::FieldContainer<double> ref_point_container( 
	ref_array_dims, const_cast<double*>(reference_point.getRawPtr()) );
    Teuchos::Array<int> phys_array_dims(3);
    phys_array_dims[0] = 1;
    phys_array_dims[1] = 1;
    phys_array_dims[2] = entity.physicalDimension();
    Intrepid::FieldContainer<double> point_container( 
	phys_array_dims, point.getRawPtr() );
    entity_cell.mapToCellPhysicalFrame( ref_point_container, point_container );
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference point.
void STKMeshEntityLocalMap::normalAtReferencePoint( 
    const Entity& entity,
    const Teuchos::ArrayView<double>& reference_point,
    const Teuchos::ArrayView<double>& normal ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//    
// Get the Intrepid cell for a given entity.
IntrepidCell& 
STKMeshEntityLocalMap::getEntityIntrepidCell( const Entity& entity ) const
{
    shards::CellTopology entity_topo = 
	stk::mesh::get_cell_topology(
	d_bulk_data->bucket( STKMeshHelpers::extractEntity(entity) ).topology() 
	    );
    DTK_CHECK( d_topo_to_cell_map.count(entity_topo.getKey()) );
    IntrepidCell& entity_cell = *d_intrepid_cells[
	d_topo_to_cell_map.find(entity_topo.getKey())->second ];
    return entity_cell;
}

//---------------------------------------------------------------------------//
void STKMeshEntityLocalMap::referenceCellCenter( 
    const Entity& entity, Intrepid::FieldContainer<double>& cell_center ) const
{
    shards::CellTopology cell_topo = 
	stk::mesh::get_cell_topology(
	d_bulk_data->bucket( STKMeshHelpers::extractEntity(entity) ).topology() 
	    );

    DTK_REQUIRE( 2 == cell_center.rank() );
    DTK_REQUIRE( Teuchos::as<unsigned>(cell_center.dimension(1)) == 
		   cell_topo.getDimension() );

    int num_cells = cell_center.dimension(0);    

    switch( cell_topo.getKey() )
    {
	case shards::Line<2>::key:
	case shards::Line<3>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = 0.0;
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
		cell_center(n,0) = 0.0;      
		cell_center(n,1) = 0.0;    
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
		cell_center(n,0) = 0.0;
		cell_center(n,1) = 0.0;
		cell_center(n,2) = 0.0;
	    }
	    break;
      
	case shards::Wedge<6>::key:
	case shards::Wedge<15>::key:
	case shards::Wedge<18>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = 1.0/3.0;
		cell_center(n,1) = 1.0/3.0;
		cell_center(n,2) = 0.0;
	    }
	    break;

	case shards::Pyramid<5>::key:
	case shards::Pyramid<13>::key:
	case shards::Pyramid<14>::key:
	    for ( int n = 0; n < num_cells; ++n )
	    {
		cell_center(n,0) = 0.0;
		cell_center(n,1) = 0.0;
		cell_center(n,2) = 1.0/4.0;
	    }
	    break;

	default:
	    DTK_INSIST( false );
	    break;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityLocalMap.cpp
//---------------------------------------------------------------------------//
