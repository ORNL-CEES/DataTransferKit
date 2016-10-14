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
 * \file tstCoarseLocalSearch.cpp
 * \author Stuart R. Slattery
 * \brief CoarseLocalSearch unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_BoxGeometry.hpp>
#include <DTK_CoarseLocalSearch.hpp>

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
TEUCHOS_UNIT_TEST( CoarseLocalSearch, coarse_local_search_test )
{
    using namespace DataTransferKit;

    // Make an entity set.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<EntitySet> entity_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );

    // Add some boxes to the set.
    int num_boxes = 5;
    for ( int i = 0; i < num_boxes; ++i )
    {
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( entity_set )
            ->addEntity( BoxGeometry( i, comm->getRank(), i, 0.0, 0.0, i, 1.0,
                                      1.0, i + 1 ) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> local_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator all_it = entity_set->entityIterator( 3 );

    // Build a coarse local search over the boxes.
    Teuchos::ParameterList plist;
    CoarseLocalSearch coarse_local_search( all_it, local_map, plist );

    // Make a point to search with.
    Teuchos::Array<double> point( 3 );
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
    plist.set<int>( "Coarse Local Search kNN", num_neighbors );
    coarse_local_search.search( point(), plist, neighbors );
    TEST_EQUALITY( num_neighbors, neighbors.size() );
    TEST_EQUALITY( 4, neighbors[0].id() );
    TEST_EQUALITY( 3, neighbors[1].id() );
}

//---------------------------------------------------------------------------//
// end tstCoarseLocalSearch.cpp
//---------------------------------------------------------------------------//
