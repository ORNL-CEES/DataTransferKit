//---------------------------------------------------------------------------//
/*! 
 * \file tstFineLocalSearch.cpp
 * \author Stuart R. Slattery
 * \brief FineLocalSearch unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_FineLocalSearch.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FineLocalSeearch, search_test_1 )
{
    using namespace DataTransferKit;

    // Add some boxes to the set.
    int num_boxes = 5;
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	boxes[i] = Box(i,0,i,0.0,0.0,i,1.0,1.0,i+1);
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

   // Build a fine local search over the boxes.
    Teuchos::ParameterList plist;
    FineLocalSearch fine_local_search( local_map );

    // Make a point to search with.
    Teuchos::Array<double> point(3);
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 2.2;

    // Search the boxes.
    Teuchos::Array<Entity> parents;
    Teuchos::Array<double> reference_coordinates;
    fine_local_search.search( boxes(), point(), plist, 
			      parents, reference_coordinates );
    TEST_EQUALITY( 1, parents.size() );
    TEST_EQUALITY( 2, parents[0].id() );
    TEST_EQUALITY( 3, reference_coordinates.size() );
    TEST_EQUALITY( point[0], reference_coordinates[0] );
    TEST_EQUALITY( point[1], reference_coordinates[1] );
    TEST_EQUALITY( point[2], reference_coordinates[2] );

    // Change the return type.
    fine_local_search.search( boxes(), point(), plist, 
			      parents, reference_coordinates );
    TEST_EQUALITY( 1, parents.size() );
    TEST_EQUALITY( 2, parents[0].id() );

    // Make a different point in no boxes.
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 5.1;
    fine_local_search.search( boxes(), point(), plist, 
			      parents, reference_coordinates );
    TEST_EQUALITY( 0, parents.size() );
    TEST_EQUALITY( 0, reference_coordinates.size() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FineLocalSeearch, search_test_2 )
{
    using namespace DataTransferKit;

    // Add some boxes to the set.
    int num_boxes = 5;
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	boxes[i] = Box(i,0,i,0.0,0.0,0.0,1.0,1.0,1.0);
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

   // Build a fine local search over the boxes.
    Teuchos::ParameterList plist;
    FineLocalSearch fine_local_search( local_map );

    // Make a point to search with.
    Teuchos::Array<double> point(3);
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 0.3;

    // Search the boxes.
    Teuchos::Array<Entity> parents;
    Teuchos::Array<double> reference_coordinates;
    fine_local_search.search( boxes(), point(), plist, 
			      parents, reference_coordinates );
    TEST_EQUALITY( 5, parents.size() );
    TEST_EQUALITY( 0, parents[0].id() );
    TEST_EQUALITY( 1, parents[1].id() );
    TEST_EQUALITY( 2, parents[2].id() );
    TEST_EQUALITY( 3, parents[3].id() );
    TEST_EQUALITY( 4, parents[4].id() );
    TEST_EQUALITY( 15, reference_coordinates.size() );
    for ( int i = 0; i < 5; ++i )
    {
	TEST_EQUALITY( point[0], reference_coordinates[3*i] );
	TEST_EQUALITY( point[1], reference_coordinates[3*i+1] );
	TEST_EQUALITY( point[2], reference_coordinates[3*i+2] );
    }

    // Make a different point in no boxes.
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 5.1;
    fine_local_search.search( boxes(), point(), plist, 
			      parents, reference_coordinates );
    TEST_EQUALITY( 0, parents.size() );
    TEST_EQUALITY( 0, reference_coordinates.size() );
}

//---------------------------------------------------------------------------//
// end tstFineLocalSearch.cpp
//---------------------------------------------------------------------------//
