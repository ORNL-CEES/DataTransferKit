//---------------------------------------------------------------------------//
/*! 
 * \file tstCoarseLocalSearch.cpp
 * \author Stuart R. Slattery
 * \brief CoarseLocalSearch unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_CoarseLocalSearch.hpp>
#include <DTK_BasicEntitySet.hpp>
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
TEUCHOS_UNIT_TEST( CoarseLocalSeearch, coarse_search_test )
{
    using namespace DataTransferKit;

    // Make an entity set.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<EntitySet> entity_set = 
	Teuchos::rcp( new BasicEntitySet(comm,3) );

    // Add some boxes to the set.
    int num_boxes = 5;
    for ( int i = 0; i < num_boxes; ++i )
    {
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(entity_set)->addEntity(
	    Box(i,comm->getRank(),i,0.0,0.0,i,1.0,1.0,i+1) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator all_it = entity_set->entityIterator( ENTITY_TYPE_VOLUME );
    
    // Build a coarse local search over the boxes.
    Teuchos::ParameterList plist;
    CoarseLocalSearch coarse_local_search( all_it, local_map, plist );

    // Make a point to search with.
    Teuchos::Array<double> point(3);
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 2.2;

    // Search the tree.
    Teuchos::Array<Entity> neighbors;
    coarse_local_search.search( point(), plist, neighbors );
    TEST_EQUALITY( num_boxes, neighbors.size() );
    TEST_EQUALITY( 2, neighbors[0].id() );
    TEST_EQUALITY( 1, neighbors[1].id() );
    TEST_EQUALITY( 3, neighbors[2].id() );
    TEST_EQUALITY( 0, neighbors[3].id() );
    TEST_EQUALITY( 4, neighbors[4].id() );

    // Make a different point.
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = 5.1;
    coarse_local_search.search( point(), plist, neighbors );
    TEST_EQUALITY( num_boxes, neighbors.size() );
    TEST_EQUALITY( 4, neighbors[0].id() );
    TEST_EQUALITY( 3, neighbors[1].id() );
    TEST_EQUALITY( 2, neighbors[2].id() );
    TEST_EQUALITY( 1, neighbors[3].id() );
    TEST_EQUALITY( 0, neighbors[4].id() );

    // Change the number of neighbors.
    int num_neighbors = 2;
    plist.set<int>("Coarse Local Search kNN", num_neighbors);
    coarse_local_search.search( point(), plist, neighbors );
    TEST_EQUALITY( num_neighbors, neighbors.size() );
    TEST_EQUALITY( 4, neighbors[0].id() );
    TEST_EQUALITY( 3, neighbors[1].id() );
}

//---------------------------------------------------------------------------//
// end tstCoarseLocalSearch.cpp
//---------------------------------------------------------------------------//
