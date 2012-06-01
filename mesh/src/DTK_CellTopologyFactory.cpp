//---------------------------------------------------------------------------//
/*!
 * \file DTK_CellTopologyFactory.cpp
 * \author Stuart R. Slattery
 * \brief Factory method defintion for cell topology data.
 */
//---------------------------------------------------------------------------//


#include "DTK_CellTopologyFactory.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_ENull.hpp>

namespace DataTransferKit
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
CellTopologyFactory::create( const moab::EntityType element_topology,
			     const int num_element_nodes )
{
    Teuchos::RCP<shards::CellTopology> new_topology;

    switch( element_topology )
    {
	case moab::MBEDGE:
	    
	    if ( num_element_nodes == 2 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Line<2> >() ) );
	    }
	    else if ( num_element_nodes == 3 )
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

	case moab::MBTRI:
	    
	    if ( num_element_nodes == 3 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Triangle<3> >() ) );
	    }
	    else if ( num_element_nodes == 4 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Triangle<4> >() ) );
	    }
	    else if ( num_element_nodes == 6 )
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

	case moab::MBQUAD:

	    if ( num_element_nodes == 4 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Quadrilateral<4> >() ) );
	    }
	    else if ( num_element_nodes == 8 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Quadrilateral<8> >() ) );
	    }
	    else if ( num_element_nodes == 9 )
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

	case moab::MBTET:
	    
	    if ( num_element_nodes == 4 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<4> >() ) );
	    }
	    else if ( num_element_nodes == 8 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<8> >() ) );
	    }
	    else if ( num_element_nodes == 10 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Tetrahedron<10> >() ) );
	    }
	    else if ( num_element_nodes == 11 )
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

	case moab::MBHEX:
	    
	    if ( num_element_nodes == 8 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Hexahedron<8> >() ) );
	    }
	    if ( num_element_nodes == 20 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Hexahedron<20> >() ) );
	    }
	    if ( num_element_nodes == 27 )
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
	    
	    testPrecondition( moab::MBEDGE == element_topology ||
			      moab::MBTRI  == element_topology ||
			      moab::MBQUAD == element_topology ||
			      moab::MBTET  == element_topology ||
			      moab::MBHEX  == element_topology ,
			      "Invalid mesh topology" );
    }

    testPostcondition( new_topology != Teuchos::null,
		       "Failure creating cell topology" );

    return new_topology;
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CellTopologyFactory.cpp
//---------------------------------------------------------------------------//

