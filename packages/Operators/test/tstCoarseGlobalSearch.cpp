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
 * \file tstCoarseGlobalSearch.cpp
 * \author Stuart R. Slattery
 * \brief CoarseGlobalSearch unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_CoarseGlobalSearch.hpp>
#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_BoxGeometry.hpp>
#include <DTK_Point.hpp>

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
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, all_to_one_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = (comm->getRank() == 0) ? 5 : 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
            BoxGeometry(i,comm_rank,i,0.0,0.0,i,1.0,1.0,i+1.0) );
    }

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    for ( int i = 0; i < num_points; ++i )
    {
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = i + 0.5;
        int id = num_points*comm_rank + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
            Point(id,comm_rank,point) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                 range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    if ( 0 == comm_rank )
    {
        int num_recv = num_points * comm_size;
        TEST_EQUALITY( num_recv, range_ids.size() );
        TEST_EQUALITY( num_recv, range_ranks.size() );
        TEST_EQUALITY( 3*num_recv, range_centroids.size() );
        for ( int i = 0; i < num_recv; ++i )
        {
            TEST_EQUALITY( range_centroids[3*i], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
            TEST_EQUALITY( range_ids[i],
                           num_points*range_ranks[i] + (range_centroids[3*i+2]-0.5) );
        }
        std::sort( range_ids.begin(), range_ids.end() );
        std::sort( range_ranks.begin(), range_ranks.end() );
        int test_id = 0;
        int idx = 0;
        for ( int i = 0; i < comm_size; ++i )
        {
            for ( int j = 0; j < num_points; ++j, ++test_id, ++idx )
            {
                TEST_EQUALITY( Teuchos::as<int>(range_ids[idx]), test_id );
                TEST_EQUALITY( range_ranks[idx], i );
            }
        }
    }
    else
    {
        TEST_EQUALITY( 0, range_ids.size() );
        TEST_EQUALITY( 0, range_ranks.size() );
        TEST_EQUALITY( 0, range_centroids.size() );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( coarse_global_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, one_to_one_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes*(comm_size-comm_rank-1) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
            BoxGeometry(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
    }

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points*comm_rank + i;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
            Point(id,comm_rank,point) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                 range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    int num_recv = num_points;
    TEST_EQUALITY( num_recv, range_ids.size() );
    TEST_EQUALITY( num_recv, range_ranks.size() );
    TEST_EQUALITY( 3*num_recv, range_centroids.size() );
    for ( int i = 0; i < num_recv; ++i )
    {
        TEST_EQUALITY( range_ranks[i], comm_size - comm_rank - 1 );
        TEST_EQUALITY( range_centroids[3*i], 0.5 );
        TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
        TEST_EQUALITY( range_centroids[3*i+2], range_ids[i]+0.5 );
    }
    std::sort( range_ids.begin(), range_ids.end() );
    for ( int j = 0; j < num_recv; ++j )
    {
        id = num_points*(comm_size-comm_rank-1) + j;
        TEST_EQUALITY( Teuchos::as<int>(range_ids[j]), id );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( coarse_global_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, no_domain_0_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();
    int comm_rank = comm->getRank();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );

    // Don't put domain entities on proc 0.
    int num_boxes = 0;
    int id = 0;
    if ( comm_rank > 0 )
    {
        num_boxes = 5;
        for ( int i = 0; i < num_boxes; ++i )
        {
            id = num_boxes*(comm_size-comm_rank-1) + i;
            Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
                BoxGeometry(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
        }
    }

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<EntityId> point_ids( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points*comm_rank + i;
        point_ids[i] = id;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
            Point(id,comm_rank,point) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                 range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    if ( comm_rank > 0 )
    {
        int num_recv = num_points;
        TEST_EQUALITY( num_recv, range_ids.size() );
        TEST_EQUALITY( num_recv, range_ranks.size() );
        TEST_EQUALITY( 3*num_recv, range_centroids.size() );
        for ( int i = 0; i < num_recv; ++i )
        {
            TEST_EQUALITY( range_ranks[i], comm_size - comm_rank - 1 );
            TEST_EQUALITY( range_centroids[3*i], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+2], range_ids[i]+0.5 );
        }
        std::sort( range_ids.begin(), range_ids.end() );
        for ( int j = 0; j < num_recv; ++j )
        {
            id = num_points*(comm_size-comm_rank-1) + j;
            TEST_EQUALITY( Teuchos::as<int>(range_ids[j]), id );
        }
    }
    else
    {
        TEST_EQUALITY( 0, range_ids.size() );
        TEST_EQUALITY( 0, range_ranks.size() );
        TEST_EQUALITY( 0, range_centroids.size() );
    }

    // Check that proc zero had all points not found.
    int num_missed = (comm_rank != comm_size-1) ? 0 : 5;
    Teuchos::Array<EntityId> missed_ids(
        coarse_global_search.getMissedRangeEntityIds() );
    TEST_EQUALITY( missed_ids.size(), num_missed );
    std::sort( point_ids.begin(), point_ids.end() );
    std::sort( missed_ids.begin(), missed_ids.end() );
    for ( int i = 0; i < num_missed; ++i )
    {
        TEST_EQUALITY( missed_ids[i], point_ids[i] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, no_range_0_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes*(comm_size-comm_rank-1) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
            BoxGeometry(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
    }

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<EntityId> point_ids( num_points );
    if ( comm_rank > 0 )
    {
        Teuchos::Array<double> point(3);
        for ( int i = 0; i < num_points; ++i )
        {
            id = num_points*comm_rank + i;
            point_ids[i] = id;
            point[0] = 0.5;
            point[1] = 0.5;
            point[2] = id + 0.5;
            Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
                Point(id,comm_rank,point) );
        }
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                 range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    if ( comm_rank < comm_size - 1 )
    {
        int num_recv = num_points;
        TEST_EQUALITY( num_recv, range_ids.size() );
        TEST_EQUALITY( num_recv, range_ranks.size() );
        TEST_EQUALITY( 3*num_recv, range_centroids.size() );
        for ( int i = 0; i < num_recv; ++i )
        {
            TEST_EQUALITY( range_ranks[i], comm_size - comm_rank - 1 );
            TEST_EQUALITY( range_centroids[3*i], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+2], range_ids[i]+0.5 );
        }
        std::sort( range_ids.begin(), range_ids.end() );
        for ( int j = 0; j < num_recv; ++j )
        {
            id = num_points*(comm_size-comm_rank-1) + j;
            TEST_EQUALITY( Teuchos::as<int>(range_ids[j]), id );
        }
    }
    else
    {
        TEST_EQUALITY( 0, range_ids.size() );
        TEST_EQUALITY( 0, range_ranks.size() );
        TEST_EQUALITY( 0, range_centroids.size() );
    }

    // Check that no missed points were found.
    int num_missed = 0;
    Teuchos::Array<EntityId> missed_ids(
        coarse_global_search.getMissedRangeEntityIds() );
    TEST_EQUALITY( missed_ids.size(), num_missed );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, many_to_many_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes*(comm_size-comm_rank-1) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
            BoxGeometry(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
    }

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 10;
    Teuchos::Array<double> point(3);
    Teuchos::Array<EntityId> point_ids( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points*comm_rank + i;
        point_ids[i] = id;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = comm_rank*5.0 + i + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
            Point(id,comm_rank,point) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                 range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    if ( comm_rank < comm_size - 1 )
    {
        int num_recv = num_points;
        TEST_EQUALITY( num_recv, range_ids.size() );
        TEST_EQUALITY( num_recv, range_ranks.size() );
        TEST_EQUALITY( 3*num_recv, range_centroids.size() );
        for ( int i = 0; i < num_recv; ++i )
        {
            TEST_EQUALITY( range_centroids[3*i], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+2],
                           range_ids[i]-5.0*range_ranks[i]+0.5 );
        }
        std::sort( range_ranks.begin(), range_ranks.end() );
        int k = 0;
        for ( int j = 1; j > -1; --j )
        {
            for ( int i = 0; i < num_points / 2; ++i, ++k )
            {
                TEST_EQUALITY( range_ranks[k], comm_size - comm_rank - 1 - j );
            }
        }
    }
    else
    {
        int num_recv = num_points / 2;
        TEST_EQUALITY( num_recv, range_ids.size() );
        TEST_EQUALITY( num_recv, range_ranks.size() );
        TEST_EQUALITY( 3*num_recv, range_centroids.size() );
        for ( int i = 0; i < num_recv; ++i )
        {
            TEST_EQUALITY( range_centroids[3*i], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
            TEST_EQUALITY( range_centroids[3*i+2],
                           range_ids[i]-5.0*range_ranks[i]+0.5 );
            TEST_EQUALITY( range_ranks[i], 0 );
        }
    }

    // Check that proc zero had some points not found.
    int num_missed = (comm_rank != comm_size-1) ? 0 : 5;
    Teuchos::Array<EntityId> missed_ids(
        coarse_global_search.getMissedRangeEntityIds() );
    TEST_EQUALITY( missed_ids.size(), num_missed );
    std::sort( point_ids.begin(), point_ids.end() );
    std::sort( missed_ids.begin(), missed_ids.end() );
    for ( int i = 0; i < num_missed; ++i )
    {
        TEST_EQUALITY( missed_ids[i], point_ids[i+5] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, point_multiple_neighbors_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int id = comm_size - comm_rank - 1;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
        BoxGeometry(id,id,id,0.0,0.0,id,1.0,1.0,id+1.0) );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    Teuchos::Array<double> point(3);
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = comm_rank;
    id = comm_rank;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
        Point(id,comm_rank,point) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                     range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    if ( comm_rank > 0 )
    {
            TEST_EQUALITY( 2, range_ids.size() );
             TEST_EQUALITY( 2, range_ranks.size() );
             TEST_EQUALITY( 6, range_centroids.size() );
             std::sort( range_ids.begin(), range_ids.end() );
             std::sort( range_ranks.begin(), range_ranks.end() );
             TEST_EQUALITY( Teuchos::as<int>(range_ids[0]), comm_size-comm_rank-1 );
             TEST_EQUALITY( Teuchos::as<int>(range_ids[1]), comm_size-comm_rank );
            TEST_EQUALITY( range_ranks[0], comm_size-comm_rank-1 );
            TEST_EQUALITY( range_ranks[1], comm_size-comm_rank );
            TEST_EQUALITY( range_centroids[0], 0.5 );
            TEST_EQUALITY( range_centroids[1], 0.5 );
            TEST_EQUALITY( range_centroids[2], comm_size-comm_rank-1 );
            TEST_EQUALITY( range_centroids[3], 0.5 );
            TEST_EQUALITY( range_centroids[4], 0.5 );
            TEST_EQUALITY( range_centroids[5], comm_size-comm_rank );
    }
    else
    {
            TEST_EQUALITY( 1, range_ids.size() );
            TEST_EQUALITY( 1, range_ranks.size() );
            TEST_EQUALITY( 3, range_centroids.size() );
            TEST_EQUALITY( Teuchos::as<int>(range_ids[0]), comm_size-comm_rank-1 );
            TEST_EQUALITY( Teuchos::as<int>(range_ranks[0]), comm_size-comm_rank-1 );
            TEST_EQUALITY( range_centroids[0], 0.5 );
            TEST_EQUALITY( range_centroids[1], 0.5 );
            TEST_EQUALITY( range_centroids[2], comm_size-comm_rank-1 );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( coarse_global_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, point_multiple_neighbors_fuzzy_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int id = comm_size - comm_rank - 1;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
        BoxGeometry(id,id,id,0.0,0.0,id,1.0,1.0,id+1.0) );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    Teuchos::Array<double> point(3);
    double fuzz = 9.9e-7;
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = comm_rank + fuzz;
    id = comm_rank;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
        Point(id,comm_rank,point) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                     range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    if ( comm_rank > 0 )
    {
            TEST_EQUALITY( 2, range_ids.size() );
             TEST_EQUALITY( 2, range_ranks.size() );
             TEST_EQUALITY( 6, range_centroids.size() );
             std::sort( range_ids.begin(), range_ids.end() );
             std::sort( range_ranks.begin(), range_ranks.end() );
             TEST_EQUALITY( Teuchos::as<int>(range_ids[0]), comm_size-comm_rank-1 );
             TEST_EQUALITY( Teuchos::as<int>(range_ids[1]), comm_size-comm_rank );
            TEST_EQUALITY( range_ranks[0], comm_size-comm_rank-1 );
            TEST_EQUALITY( range_ranks[1], comm_size-comm_rank );
            TEST_EQUALITY( range_centroids[0], 0.5 );
            TEST_EQUALITY( range_centroids[1], 0.5 );
            TEST_EQUALITY( range_centroids[2], comm_size-comm_rank-1 + fuzz);
            TEST_EQUALITY( range_centroids[3], 0.5 );
            TEST_EQUALITY( range_centroids[4], 0.5 );
            TEST_EQUALITY( range_centroids[5], comm_size-comm_rank + fuzz);
    }
    else
    {
            TEST_EQUALITY( 1, range_ids.size() );
            TEST_EQUALITY( 1, range_ranks.size() );
            TEST_EQUALITY( 3, range_centroids.size() );
            TEST_EQUALITY( Teuchos::as<int>(range_ids[0]), comm_size-comm_rank-1 );
            TEST_EQUALITY( Teuchos::as<int>(range_ranks[0]), comm_size-comm_rank-1 );
            TEST_EQUALITY( range_centroids[0], 0.5 );
            TEST_EQUALITY( range_centroids[1], 0.5 );
            TEST_EQUALITY( range_centroids[2], comm_size-comm_rank-1 + fuzz);
    }

    // Check that no missed points were found.
    TEST_EQUALITY( coarse_global_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, point_single_neighbor_fuzzy_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int id = comm_size - comm_rank - 1;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
        BoxGeometry(id,id,id,0.0,0.0,id,1.0,1.0,id+1.0) );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    Teuchos::Array<double> point(3);
    double fuzz = 2.0e-6;
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = comm_rank - fuzz;
    id = comm_rank;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
        Point(id,comm_rank,point) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                     range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    if ( comm_rank > 0 )
    {
            TEST_EQUALITY( 1, range_ids.size() );
            TEST_EQUALITY( 1, range_ranks.size() );
            TEST_EQUALITY( 3, range_centroids.size() );
            TEST_EQUALITY( Teuchos::as<int>(range_ids[0]), comm_size-comm_rank );
            TEST_EQUALITY( Teuchos::as<int>(range_ranks[0]), comm_size-comm_rank );
            TEST_EQUALITY( range_centroids[0], 0.5 );
            TEST_EQUALITY( range_centroids[1], 0.5 );
            TEST_EQUALITY( range_centroids[2], comm_size-comm_rank - fuzz);
        TEST_EQUALITY( coarse_global_search.getMissedRangeEntityIds().size(), 0 );
    }
    else
    {
        TEST_EQUALITY( 0, range_ids.size() );
            TEST_EQUALITY( 0, range_ranks.size() );
            TEST_EQUALITY( 0, range_centroids.size() );
        TEST_EQUALITY( coarse_global_search.getMissedRangeEntityIds().size(), 1 );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CoarseGlobalSearch, missed_range_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes*(comm_size-comm_rank-1) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
            BoxGeometry(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
    }

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a coarse global search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>("Track Missed Range Entities",true);
    CoarseGlobalSearch coarse_global_search( comm, 3, domain_it, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points*comm_rank + i;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
            Point(id,comm_rank,point) );
    }

    // Add a bad point.
    id = num_points*comm_rank + 1000;
    point[0] = -100.0;
    point[1] = 0.0;
    point[2] = 0.0;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
        Point(id,comm_rank,point) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Search the tree.
    Teuchos::Array<EntityId> range_ids;
    Teuchos::Array<int> range_ranks;
    Teuchos::Array<double> range_centroids;
    coarse_global_search.search( range_it, range_map, plist,
                                 range_ids, range_ranks, range_centroids );

    // Check the results of the search.
    int num_recv = num_points;
    TEST_EQUALITY( num_recv, range_ids.size() );
    TEST_EQUALITY( num_recv, range_ranks.size() );
    TEST_EQUALITY( 3*num_recv, range_centroids.size() );
    for ( int i = 0; i < num_recv; ++i )
    {
        TEST_EQUALITY( range_ranks[i], comm_size - comm_rank - 1 );
        TEST_EQUALITY( range_centroids[3*i], 0.5 );
        TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
        TEST_EQUALITY( range_centroids[3*i+2], range_ids[i]+0.5 );
    }
    std::sort( range_ids.begin(), range_ids.end() );
    for ( int j = 0; j < num_recv; ++j )
    {
        id = num_points*(comm_size-comm_rank-1) + j;
        TEST_EQUALITY( Teuchos::as<int>(range_ids[j]), id );
    }

    // Check that the bad point was found.
    Teuchos::ArrayView<const EntityId> missed_range =
        coarse_global_search.getMissedRangeEntityIds();
    TEST_EQUALITY( missed_range.size(), 1 );
    TEST_EQUALITY( missed_range[0],
                   Teuchos::as<EntityId>(num_points*comm_rank + 1000) );
}

//---------------------------------------------------------------------------//
// end tstCoarseGlobalSearch.cpp
//---------------------------------------------------------------------------//
