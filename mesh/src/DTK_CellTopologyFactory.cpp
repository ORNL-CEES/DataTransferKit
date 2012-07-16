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
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
CellTopologyFactory::CellTopologyFactory()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
CellTopologyFactory::~CellTopologyFactory()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Factory method. If the number of element nodes provided is not
 * supported for that topology, the topology for the linear element will be
 * created. 
 */
CellTopologyFactory::RCP_CellTopology
CellTopologyFactory::create( const moab::EntityType element_topology,
			     const int num_element_nodes )
{
    RCP_CellTopology new_topology;

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

	case moab::MBPYRAMID:
	    
	    if ( num_element_nodes == 5 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Pyramid<5> >() ) );
	    }
	    if ( num_element_nodes == 13 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Pyramid<13> >() ) );
	    }
	    if ( num_element_nodes == 14 )
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Pyramid<14> >() ) );
	    }
	    else 
	    {
		new_topology = Teuchos::rcp( 
		    new shards::CellTopology(
			shards::getCellTopologyData< shards::Pyramid<> >() ) );
	    }
	    break;

	default:
	    
	    testPrecondition( moab::MBEDGE    == element_topology ||
			      moab::MBTRI     == element_topology ||
			      moab::MBQUAD    == element_topology ||
			      moab::MBTET     == element_topology ||
			      moab::MBHEX     == element_topology ||
			      moab::MBPYRAMID == element_topology,
			      "Invalid mesh topology" );
    }

    testPostcondition( new_topology != Teuchos::null,
		       "Failure creating cell topology" );

    return new_topology;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CellTopologyFactory.cpp
//---------------------------------------------------------------------------//

