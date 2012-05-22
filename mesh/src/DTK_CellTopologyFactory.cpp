//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file CellTopologyFactory.cpp
 * \author Stuart Slattery
 * \brief Factory method defintion for cell topology data.
 */
//---------------------------------------------------------------------------//

#include "Exception.hpp"
#include "CellTopologyFactory.hpp"

#include <iMesh.h>

#include <Teuchos_ENull.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
CellTopologyFactory::CellTopologyFactory()
{ /* ... */ }

/*!
 * \brief Destructor.
 */
CellTopologyFactory::~CellTopologyFactory()
{ /* ... */ }

/*!
 * \brief Factory method.
 */
Teuchos::RCP<shards::CellTopology> 
CellTopologyFactory::create( const int entity_topology,
			     const int num_entity_nodes )
{
    Teuchos::RCP<shards::CellTopology> new_topology;

    switch( entity_topology )
    {
	case iMesh_LINE_SEGMENT:
	    
	    if ( num_entity_nodes == 2 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Line<2> >() ) );
	    }
	    else if ( num_entity_nodes == 3 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Line<3> >() ) );
	    }
	    else {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Line<> >() ) );
	    }
	    break;

	case iMesh_TRIANGLE:
	    
	    if ( num_entity_nodes == 3 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Triangle<3> >() ) );
	    }
	    else if ( num_entity_nodes == 4 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Triangle<4> >() ) );
	    }
	    else if ( num_entity_nodes == 6 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Triangle<6> >() ) );
	    }
	    else
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Triangle<> >() ) );
	    }
	    break;

	case iMesh_QUADRILATERAL:

	    if ( num_entity_nodes == 4 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Quadrilateral<4> >() ) );
	    }
	    else if ( num_entity_nodes == 8 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Quadrilateral<8> >() ) );
	    }
	    else if ( num_entity_nodes == 9 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Quadrilateral<9> >() ) );
	    }
	    else
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Quadrilateral<> >() ) );
	    }
	    break;

	case iMesh_TETRAHEDRON:
	    
	    if ( num_entity_nodes == 4 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<4> >() ) );
	    }
	    else if ( num_entity_nodes == 8 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<8> >() ) );
	    }
	    else if ( num_entity_nodes == 10 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<10> >() ) );
	    }
	    else if ( num_entity_nodes == 11 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<11> >() ) );
	    }
	    else 
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<> >() ) );
	    }
	    break;

	case iMesh_HEXAHEDRON:
	    
	    if ( num_entity_nodes == 8 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Hexahedron<8> >() ) );
	    }
	    if ( num_entity_nodes == 20 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Hexahedron<20> >() ) );
	    }
	    if ( num_entity_nodes == 27 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Hexahedron<27> >() ) );
	    }
	    else 
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Hexahedron<> >() ) );
	    }
	    break;

	default:
	    
	    testPrecondition( iMesh_LINE_SEGMENT  == entity_topology ||
			      iMesh_TRIANGLE      == entity_topology ||
			      iMesh_QUADRILATERAL == entity_topology ||
			      iMesh_TETRAHEDRON   == entity_topology ||
			      iMesh_HEXAHEDRON    == entity_topology ,
			      "Invalid mesh topology" );
    }

    testPostcondition( new_topology != Teuchos::null,
		       "Failure creating cell topology" );

    return new_topology;
}

}

//---------------------------------------------------------------------------//
// end CellTopologyFactory.cpp
//---------------------------------------------------------------------------//

