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
 * \file   tstReferenceHexShapeFunction.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Nodal shape function test.
 */
//---------------------------------------------------------------------------//

#include "reference_implementation/DTK_ReferenceHexShapeFunction.hpp"
#include "reference_implementation/DTK_ReferenceHex.hpp"
#include "reference_implementation/DTK_ReferenceNode.hpp"

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>

//---------------------------------------------------------------------------//
// Hex test.
TEUCHOS_UNIT_TEST( ReferenceHexShapeFunction, hex_test )
{
    // Create the nodes;
    int num_nodes = 8;
    Teuchos::Array<DataTransferKit::Entity> nodes( num_nodes );
    nodes[0] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 0, 0.0, 0.0, 0.0 );
    nodes[1] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 1, 2.0, 0.0, 0.0 );
    nodes[2] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 2, 2.0, 2.0, 0.0 );
    nodes[3] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 3, 0.0, 2.0, 0.0 );
    nodes[4] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 4, 0.0, 0.0, 2.0 );
    nodes[5] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 5, 2.0, 0.0, 2.0 );
    nodes[6] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 6, 2.0, 2.0, 2.0 );
    nodes[7] =
        DataTransferKit::UnitTest::ReferenceNode( 0, 7, 0.0, 2.0, 2.0 );
    
    // Make a hex.
    DataTransferKit::Entity hex = 
        DataTransferKit::UnitTest::ReferenceHex( 0, 0, nodes );
    
    // Create a shape function.
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> shape_function =
	Teuchos::rcp( new DataTransferKit::UnitTest::ReferenceHexShapeFunction() );

    // Test the shape function dof ids for the hex.
    Teuchos::Array<DataTransferKit::SupportId> dof_ids;
    shape_function->entitySupportIds( hex, dof_ids );
    TEST_EQUALITY( num_nodes, dof_ids.size() );
    for ( int n = 0; n < num_nodes; ++n )
    {
	TEST_EQUALITY( dof_ids[n], nodes[n].id() );
    }

    // Test the value evaluation for the hex.
    Teuchos::Array<double> ref_point( 3, 0.0 );
    Teuchos::Array<double> values;
    shape_function->evaluateValue( hex, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    for ( int n = 0; n < num_nodes; ++n )
    {
	TEST_EQUALITY( values[n], 1.0 / num_nodes );
    }
    ref_point[0] = -1.0;
    ref_point[1] = -1.0;
    ref_point[2] = -1.0;
    shape_function->evaluateValue( hex, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    TEST_EQUALITY( values[0], 1.0 );
    for ( int n = 1; n < num_nodes; ++n )
    {
	TEST_EQUALITY( values[n], 0.0 );
    }    

    // Test the gradient evaluation for the hex.
    Teuchos::Array<Teuchos::Array<double> > grads;
    ref_point.assign( 3, 0.0 );
    shape_function->evaluateGradient( hex, ref_point(), grads );
    TEST_EQUALITY( grads.size(), num_nodes );
    for ( int n = 0; n < num_nodes; ++n )
    {
	TEST_EQUALITY( Teuchos::as<int>(grads[n].size()), 3 );
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

    // Test the shape function dof ids for the nodes.
    for ( int n = 0; n < num_nodes; ++n )
    {
	dof_ids.clear();
	shape_function->entitySupportIds( nodes[n], dof_ids );
	TEST_EQUALITY( dof_ids.size(), 1 );
	TEST_EQUALITY( dof_ids[0], nodes[n].id() );
    }
}

//---------------------------------------------------------------------------//
// end of tstReferenceHexShapeFunction.cpp
//---------------------------------------------------------------------------//
