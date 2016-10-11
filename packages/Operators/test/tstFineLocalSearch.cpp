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
 * \file tstFineLocalSearch.cpp
 * \author Stuart R. Slattery
 * \brief FineLocalSearch unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_BoxGeometry.hpp>
#include <DTK_FineLocalSearch.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FineLocalSearch, search_test_1 )
{
    using namespace DataTransferKit;

    // Add some boxes to the set.
    int num_boxes = 5;
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
        boxes[i] = BoxGeometry( i, 0, i, 0.0, 0.0, i, 1.0, 1.0, i + 1 );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> local_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Build a fine local search over the boxes.
    Teuchos::ParameterList plist;
    FineLocalSearch fine_local_search( local_map );

    // Make a point to search with.
    Teuchos::Array<double> point( 3 );
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 2.2;

    // Search the boxes.
    Teuchos::Array<Entity> parents;
    Teuchos::Array<double> reference_coordinates;
    fine_local_search.search( boxes(), point(), plist, parents,
                              reference_coordinates );
    TEST_EQUALITY( 1, parents.size() );
    TEST_EQUALITY( 2, parents[0].id() );
    TEST_EQUALITY( 3, reference_coordinates.size() );
    TEST_EQUALITY( point[0], reference_coordinates[0] );
    TEST_EQUALITY( point[1], reference_coordinates[1] );
    TEST_EQUALITY( point[2], reference_coordinates[2] );

    // Change the return type.
    fine_local_search.search( boxes(), point(), plist, parents,
                              reference_coordinates );
    TEST_EQUALITY( 1, parents.size() );
    TEST_EQUALITY( 2, parents[0].id() );

    // Make a different point in no boxes.
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 5.1;
    fine_local_search.search( boxes(), point(), plist, parents,
                              reference_coordinates );
    TEST_EQUALITY( 0, parents.size() );
    TEST_EQUALITY( 0, reference_coordinates.size() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FineLocalSearch, search_test_2 )
{
    using namespace DataTransferKit;

    // Add some boxes to the set.
    int num_boxes = 5;
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
        boxes[i] = BoxGeometry( i, 0, i, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> local_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Build a fine local search over the boxes.
    Teuchos::ParameterList plist;
    FineLocalSearch fine_local_search( local_map );

    // Make a point to search with.
    Teuchos::Array<double> point( 3 );
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 0.3;

    // Search the boxes.
    Teuchos::Array<Entity> parents;
    Teuchos::Array<double> reference_coordinates;
    fine_local_search.search( boxes(), point(), plist, parents,
                              reference_coordinates );
    TEST_EQUALITY( 5, parents.size() );
    TEST_EQUALITY( 0, parents[0].id() );
    TEST_EQUALITY( 1, parents[1].id() );
    TEST_EQUALITY( 2, parents[2].id() );
    TEST_EQUALITY( 3, parents[3].id() );
    TEST_EQUALITY( 4, parents[4].id() );
    TEST_EQUALITY( 15, reference_coordinates.size() );
    for ( int i = 0; i < 5; ++i )
    {
        TEST_EQUALITY( point[0], reference_coordinates[3 * i] );
        TEST_EQUALITY( point[1], reference_coordinates[3 * i + 1] );
        TEST_EQUALITY( point[2], reference_coordinates[3 * i + 2] );
    }

    // Make a different point in no boxes.
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 5.1;
    fine_local_search.search( boxes(), point(), plist, parents,
                              reference_coordinates );
    TEST_EQUALITY( 0, parents.size() );
    TEST_EQUALITY( 0, reference_coordinates.size() );
}

//---------------------------------------------------------------------------//
// end tstFineLocalSearch.cpp
//---------------------------------------------------------------------------//
