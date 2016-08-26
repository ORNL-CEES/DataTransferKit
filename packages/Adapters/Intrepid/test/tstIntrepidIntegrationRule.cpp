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
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstIntrepidIntegrationRule.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Integration rule function test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_IntrepidIntegrationRule.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <Teuchos_DefaultMpiComm.hpp>

#include <Shards_BasicTopologies.hpp>

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( IntrepidIntegrationRule, hex_8_test )
{
    // Create an integration rule.
    DataTransferKit::IntrepidIntegrationRule integration_rule;

    // Create a cell topology.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();    
    
    // Test the integration rule.
    Teuchos::Array<Teuchos::Array<double> > p_1;
    Teuchos::Array<double> w_1;
    integration_rule.getIntegrationRule( element_topo, 1, p_1, w_1 );
    TEST_EQUALITY( 1, w_1.size() );
    TEST_EQUALITY( 1, p_1.size() );
    TEST_EQUALITY( 3, p_1[0].size() );
    TEST_EQUALITY( 8.0, w_1[0] );
    TEST_EQUALITY( 0.0, p_1[0][0] );
    TEST_EQUALITY( 0.0, p_1[0][1] );
    TEST_EQUALITY( 0.0, p_1[0][2] );

    Teuchos::Array<Teuchos::Array<double> > p_2;
    Teuchos::Array<double> w_2;
    integration_rule.getIntegrationRule( element_topo, 2, p_2, w_2 );
    TEST_EQUALITY( 8, w_2.size() );
    TEST_EQUALITY( 8, p_2.size() );
    for ( int i = 0; i < 8; ++i )
    {
	TEST_EQUALITY( w_2[i], 1.0 );
	TEST_EQUALITY( p_2[i].size(), 3 );

	for ( int d = 0; d < 3; ++d )
	{
	    TEST_FLOATING_EQUALITY(
		std::abs(p_2[i][d]), 1.0 / std::sqrt(3.0), 1.0e-15 );
	}
    }
}

//---------------------------------------------------------------------------//
// end of tstIntrepidIntegrationRule.cpp
//---------------------------------------------------------------------------//
