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
 * \file tstReferenceHexLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief ReferenceHexLocalMap unit tests.
 */
//---------------------------------------------------------------------------//

#include "reference_implementation/DTK_ReferenceHexLocalMap.hpp"
#include "reference_implementation/DTK_ReferenceHex.hpp"
#include "reference_implementation/DTK_ReferenceNode.hpp"

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ReferenceHexLocalMap, local_map_test )
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
    
    // Create a local map from the bulk data.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
	Teuchos::rcp( new DataTransferKit::UnitTest::ReferenceHexLocalMap() );

    // Test the measure.
    TEST_EQUALITY( local_map->measure(hex), 8.0 );
    TEST_EQUALITY( local_map->measure(nodes[0]), 0.0 );

    // Test the centroid.
    Teuchos::Array<double> centroid( 3, 0.0 );
    local_map->centroid( hex, centroid() );
    TEST_EQUALITY( centroid[0], 1.0 );
    TEST_EQUALITY( centroid[1], 1.0 );
    TEST_EQUALITY( centroid[2], 1.0 );
    local_map->centroid( nodes[0], centroid() );
    TEST_EQUALITY( centroid[0], 0.0 );
    TEST_EQUALITY( centroid[1], 0.0 );
    TEST_EQUALITY( centroid[2], 0.0 );    
    local_map->centroid( nodes[1], centroid() );
    TEST_EQUALITY( centroid[0], 2.0 );
    TEST_EQUALITY( centroid[1], 0.0 );
    TEST_EQUALITY( centroid[2], 0.0 );    
    local_map->centroid( nodes[2], centroid() );
    TEST_EQUALITY( centroid[0], 2.0 );
    TEST_EQUALITY( centroid[1], 2.0 );
    TEST_EQUALITY( centroid[2], 0.0 );
    local_map->centroid( nodes[3], centroid() );
    TEST_EQUALITY( centroid[0], 0.0 );
    TEST_EQUALITY( centroid[1], 2.0 );
    TEST_EQUALITY( centroid[2], 0.0 );        
    local_map->centroid( nodes[4], centroid() );
    TEST_EQUALITY( centroid[0], 0.0 );
    TEST_EQUALITY( centroid[1], 0.0 );
    TEST_EQUALITY( centroid[2], 2.0 );    
    local_map->centroid( nodes[5], centroid() );
    TEST_EQUALITY( centroid[0], 2.0 );
    TEST_EQUALITY( centroid[1], 0.0 );
    TEST_EQUALITY( centroid[2], 2.0 );    
    local_map->centroid( nodes[6], centroid() );
    TEST_EQUALITY( centroid[0], 2.0 );
    TEST_EQUALITY( centroid[1], 2.0 );
    TEST_EQUALITY( centroid[2], 2.0 );    
    local_map->centroid( nodes[7], centroid() );
    TEST_EQUALITY( centroid[0], 0.0 );
    TEST_EQUALITY( centroid[1], 2.0 );
    TEST_EQUALITY( centroid[2], 2.0 );    
    
    // Make a good point and a bad point.
    Teuchos::Array<double> good_point( 3 );
    good_point[0] = 0.5;
    good_point[1] = 1.5;
    good_point[2] = 1.0;
    Teuchos::Array<double> bad_point( 3 );
    bad_point[0] = 0.75;
    bad_point[1] = -1.75;
    bad_point[2] = 0.35;

    // Test the reference frame safeguard.
    TEST_ASSERT(
    	local_map->isSafeToMapToReferenceFrame(hex,good_point()) );
    TEST_ASSERT(
    	!local_map->isSafeToMapToReferenceFrame(hex,bad_point()) );

    // Test the mapping to reference frame.
    Teuchos::Array<double> ref_good_point( 3 );
    bool good_map = local_map->mapToReferenceFrame( 
    	hex, good_point(), ref_good_point() );
    TEST_ASSERT( good_map );
    TEST_EQUALITY( ref_good_point[0], -0.5 );
    TEST_EQUALITY( ref_good_point[1], 0.5 );
    TEST_EQUALITY( ref_good_point[2], 0.0 );
			    
    Teuchos::Array<double> ref_bad_point( 3 );
    bool bad_map = local_map->mapToReferenceFrame( 
    	hex, bad_point(), ref_bad_point() );
    TEST_ASSERT( bad_map );

    // Test the point inclusion.
    TEST_ASSERT( local_map->checkPointInclusion(hex,ref_good_point()) );
    TEST_ASSERT( !local_map->checkPointInclusion(hex,ref_bad_point()) );

    // Test the map to physical frame.
    Teuchos::Array<double> phy_good_point( 3 );
    local_map->mapToPhysicalFrame(hex,ref_good_point(),phy_good_point());
    TEST_EQUALITY( good_point[0], phy_good_point[0] );
    TEST_EQUALITY( good_point[1], phy_good_point[1] );
    TEST_EQUALITY( good_point[2], phy_good_point[2] );

    Teuchos::Array<double> phy_bad_point( 3 );
    local_map->mapToPhysicalFrame(hex,ref_bad_point(),phy_bad_point());
    TEST_EQUALITY( bad_point[0], phy_bad_point[0] );
    TEST_EQUALITY( bad_point[1], phy_bad_point[1] );
    TEST_EQUALITY( bad_point[2], phy_bad_point[2] );
}

//---------------------------------------------------------------------------//
// end tstReferenceHexLocalMap.cpp
//---------------------------------------------------------------------------//

