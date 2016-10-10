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
 * \file   tstIntrepidShapeFunction.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Nodal shape function test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_IntrepidShapeFunction.hpp"

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
TEUCHOS_UNIT_TEST( IntrepidShapeFunction, hex_8_test )
{
    // Create a shape function.
    DataTransferKit::IntrepidShapeFunction shape_function;

    // Create a cell topology.
    shards::CellTopology element_topo =
        shards::getCellTopologyData<shards::Hexahedron<8> >();

    // Test the value evaluation for the hex.
    int space_dim = 3;
    int num_nodes = 8;
    Teuchos::Array<double> ref_point( space_dim, 0.0 );
    Teuchos::Array<double> values;
    shape_function.evaluateValue( element_topo, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    for ( int n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( values[n], 1.0 / num_nodes );
    }
    ref_point[0] = -1.0;
    ref_point[1] = -1.0;
    ref_point[2] = -1.0;
    shape_function.evaluateValue( element_topo, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    TEST_EQUALITY( values[0], 1.0 );
    for ( int n = 1; n < num_nodes; ++n )
    {
        TEST_EQUALITY( values[n], 0.0 );
    }

    // Test the gradient evaluation for the hex.
    Teuchos::Array<Teuchos::Array<double> > grads;
    ref_point.assign( 3, 0.0 );
    shape_function.evaluateGradient( element_topo, ref_point(), grads );
    TEST_EQUALITY( grads.size(), num_nodes );
    for ( int n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( Teuchos::as<int>(grads[n].size()), space_dim );
    }

    TEST_EQUALITY( grads[0][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[0][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[0][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[1][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[1][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[1][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[2][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[2][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[2][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[3][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[3][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[3][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[4][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[4][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[4][2], 1.0 / num_nodes );

    TEST_EQUALITY( grads[5][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[5][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[5][2], 1.0 / num_nodes );

    TEST_EQUALITY( grads[6][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[6][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[6][2], 1.0 / num_nodes );

    TEST_EQUALITY( grads[7][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[7][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[7][2], 1.0 / num_nodes );
}

//---------------------------------------------------------------------------//
// end of tstIntrepidShapeFunction.cpp
//---------------------------------------------------------------------------//
