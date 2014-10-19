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
    Teuchos::Array<Entity> mapped_boxes;
    fine_local_search.search( boxes(), point(), plist, mapped_boxes );
    TEST_EQUALITY( 1, mapped_boxes.size() );
    TEST_EQUALITY( 2, mapped_boxes[0].id() );

    // Change the return type.
    plist.set<bool>("Fine Local Search Return All", true);
    fine_local_search.search( boxes(), point(), plist, mapped_boxes );
    TEST_EQUALITY( 1, mapped_boxes.size() );
    TEST_EQUALITY( 2, mapped_boxes[0].id() );

    // Make a different point in no boxes.
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 5.1;
    fine_local_search.search( boxes(), point(), plist, mapped_boxes );
    TEST_EQUALITY( 0, mapped_boxes.size() );
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
    Teuchos::Array<Entity> mapped_boxes;
    fine_local_search.search( boxes(), point(), plist, mapped_boxes );
    TEST_EQUALITY( 1, mapped_boxes.size() );
    TEST_EQUALITY( 0, mapped_boxes[0].id() );

    // Change the return type.
    plist.set<bool>("Fine Local Search Return All", true);
    fine_local_search.search( boxes(), point(), plist, mapped_boxes );
    TEST_EQUALITY( 5, mapped_boxes.size() );
    TEST_EQUALITY( 0, mapped_boxes[0].id() );
    TEST_EQUALITY( 1, mapped_boxes[1].id() );
    TEST_EQUALITY( 2, mapped_boxes[2].id() );
    TEST_EQUALITY( 3, mapped_boxes[3].id() );
    TEST_EQUALITY( 4, mapped_boxes[4].id() );

    // Make a different point in no boxes.
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 5.1;
    fine_local_search.search( boxes(), point(), plist, mapped_boxes );
    TEST_EQUALITY( 0, mapped_boxes.size() );
}

//---------------------------------------------------------------------------//
// end tstFineLocalSearch.cpp
//---------------------------------------------------------------------------//
