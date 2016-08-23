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
 * \file tstReferenceHex.cpp
 * \author Stuart R. Slattery
 * \brief ReferenceHex unit tests.
 */
//---------------------------------------------------------------------------//

#include "reference_implementation/DTK_ReferenceHex.hpp"
#include "reference_implementation/DTK_ReferenceNode.hpp"

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <iostream>

//---------------------------------------------------------------------------//
// Hex test.
TEUCHOS_UNIT_TEST( ReferenceHex, hex_test )
{
    // Get the comm.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    
    // Create the nodes;
    int num_nodes = 8;
    Teuchos::Array<DataTransferKit::Entity> nodes( num_nodes );
    nodes[0] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 0, 0.0, 0.0, 0.0 );
    nodes[1] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 1, 1.0, 0.0, 0.0 );
    nodes[2] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 2, 1.0, 1.0, 0.0 );
    nodes[3] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 3, 0.0, 1.0, 0.0 );
    nodes[4] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 4, 0.0, 0.0, 1.0 );
    nodes[5] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 5, 1.0, 0.0, 1.0 );
    nodes[6] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 6, 1.0, 1.0, 1.0 );
    nodes[7] =
        DataTransferKit::UnitTest::ReferenceNode( comm_rank, 7, 0.0, 1.0, 1.0 );
    
    // Make a hex.
    int hex_id = 23 + comm_rank;
    DataTransferKit::Entity entity = 
        DataTransferKit::UnitTest::ReferenceHex( hex_id, comm_rank, nodes );

    // Print out the entity.
    std::cout << entity.description() << std::endl;
    Teuchos::RCP<Teuchos::FancyOStream>
	fancy_out = Teuchos::VerboseObjectBase::getDefaultOStream();
    entity.describe( *fancy_out );
    
    // Test the entity.
    TEST_EQUALITY( hex_id, entity.id() );
    TEST_EQUALITY( comm_rank, entity.ownerRank() );
    TEST_EQUALITY( 3, entity.topologicalDimension() );
    TEST_EQUALITY( 3, entity.physicalDimension() );
    TEST_ASSERT( !entity.inBlock(0) );
    TEST_ASSERT( !entity.onBoundary(0) );

    Teuchos::Tuple<double,6> hex_bounds;
    entity.boundingBox( hex_bounds );
    TEST_EQUALITY( 0.0, hex_bounds[0] );
    TEST_EQUALITY( 0.0, hex_bounds[1] );
    TEST_EQUALITY( 0.0, hex_bounds[2] );
    TEST_EQUALITY( 1.0, hex_bounds[3] );
    TEST_EQUALITY( 1.0, hex_bounds[4] );
    TEST_EQUALITY( 1.0, hex_bounds[5] );
}

//---------------------------------------------------------------------------//
// end tstReferenceHex.cpp
//---------------------------------------------------------------------------//

