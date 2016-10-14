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
 * \file tstParallelSearch.cpp
 * \author Stuart R. Slattery
 * \brief ParallelSearch unit tests.
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
#include <DTK_ParallelSearch.hpp>
#include <DTK_Point.hpp>

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
TEUCHOS_UNIT_TEST( ParallelSearch, all_to_one_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_boxes = ( comm->getRank() == 0 ) ? 5 : 0;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
        box_ids[i] = i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
            ->addEntity( BoxGeometry( i, comm_rank, i, 0.0, 0.0, i, 1.0, 1.0,
                                      i + 1.0 ) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_points = 5;
    Teuchos::Array<double> point( 3 );
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = i + 0.5;
        point_ids[i] = num_points * comm_rank + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
            ->addEntity( Point( point_ids[i], comm_rank, point ) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    Teuchos::Array<EntityId> local_range;
    Teuchos::Array<EntityId> local_domain;
    Teuchos::Array<EntityId> range_entities;
    for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
          ++domain_it )
    {
        parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                    range_entities );
        TEST_EQUALITY( comm_size, range_entities.size() );
        for ( int n = 0; n < comm_size; ++n )
        {
            TEST_EQUALITY( range_entities[n] % num_points, domain_it->id() );
            local_range.push_back( range_entities[n] );
            local_domain.push_back( domain_it->id() );
        }
    }

    int range_size = ( comm_rank == 0 ) ? num_points * comm_size : 0;
    TEST_EQUALITY( local_range.size(), range_size );

    Teuchos::ArrayView<const double> range_coords;
    Teuchos::Array<EntityId> domain_entities;
    for ( int i = 0; i < num_points; ++i )
    {
        parallel_search.getDomainEntitiesFromRange( point_ids[i],
                                                    domain_entities );
        TEST_EQUALITY( 1, domain_entities.size() );
        TEST_EQUALITY( point_ids[i] % num_points, domain_entities[0] );
    }

    for ( int i = 0; i < range_size; ++i )
    {
        parallel_search.rangeParametricCoordinatesInDomain(
            local_domain[i], local_range[i], range_coords );
        TEST_EQUALITY( range_coords[0], 0.5 );
        TEST_EQUALITY( range_coords[1], 0.5 );
        TEST_EQUALITY( range_coords[2], local_domain[i] + 0.5 );
        TEST_EQUALITY( Teuchos::as<int>( ( local_range[i] - local_domain[i] ) /
                                         num_points ),
                       parallel_search.rangeEntityOwnerRank( local_range[i] ) );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( parallel_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, one_to_one_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes * ( comm_size - comm_rank - 1 ) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
            ->addEntity( BoxGeometry( id, comm_rank, id, 0.0, 0.0, id, 1.0, 1.0,
                                      id + 1.0 ) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_points = 5;
    Teuchos::Array<double> point( 3 );
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points * comm_rank + i;
        point_ids[i] = id;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
            ->addEntity( Point( id, comm_rank, point ) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    Teuchos::Array<EntityId> local_range;
    Teuchos::Array<EntityId> local_domain;
    Teuchos::Array<EntityId> range_entities;
    for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
          ++domain_it )
    {
        parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                    range_entities );
        TEST_EQUALITY( 1, range_entities.size() );
        TEST_EQUALITY( range_entities[0], domain_it->id() );
        local_range.push_back( range_entities[0] );
        local_domain.push_back( domain_it->id() );
    }
    TEST_EQUALITY( local_range.size(), 5 );
    Teuchos::ArrayView<const double> range_coords;
    Teuchos::Array<EntityId> domain_entities;
    for ( int i = 0; i < 5; ++i )
    {
        parallel_search.getDomainEntitiesFromRange( point_ids[i],
                                                    domain_entities );
        TEST_EQUALITY( 1, domain_entities.size() );
        TEST_EQUALITY( point_ids[i], domain_entities[0] );
    }
    for ( int i = 0; i < 5; ++i )
    {
        parallel_search.rangeParametricCoordinatesInDomain(
            local_domain[i], local_range[i], range_coords );
        TEST_EQUALITY( range_coords[0], 0.5 );
        TEST_EQUALITY( range_coords[1], 0.5 );
        TEST_EQUALITY( range_coords[2], local_range[i] + 0.5 );
        TEST_EQUALITY( comm_size - comm_rank - 1,
                       parallel_search.rangeEntityOwnerRank( local_range[i] ) );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( parallel_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, no_domain_0_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();
    int comm_rank = comm->getRank();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );

    // Don't put domain entities on proc 0.
    int num_boxes = 0;
    int id = 0;
    if ( comm_rank > 0 )
    {
        num_boxes = 5;
        for ( int i = 0; i < num_boxes; ++i )
        {
            id = num_boxes * ( comm_size - comm_rank - 1 ) + i;
            Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
                ->addEntity( BoxGeometry( id, comm_rank, id, 0.0, 0.0, id, 1.0,
                                          1.0, id + 1.0 ) );
        }
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_points = 5;
    Teuchos::Array<double> point( 3 );
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points * comm_rank + i;
        point_ids[i] = id;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
            ->addEntity( Point( id, comm_rank, point ) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    if ( comm_rank > 0 )
    {
        Teuchos::Array<EntityId> local_range;
        Teuchos::Array<EntityId> local_domain;
        Teuchos::Array<EntityId> range_entities;
        for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
              ++domain_it )
        {
            parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                        range_entities );
            TEST_EQUALITY( 1, range_entities.size() );
            TEST_EQUALITY( range_entities[0], domain_it->id() );
            local_range.push_back( range_entities[0] );
            local_domain.push_back( domain_it->id() );
        }
        TEST_EQUALITY( local_range.size(), 5 );
        Teuchos::ArrayView<const double> range_coords;
        Teuchos::Array<EntityId> domain_entities;
        if ( comm_rank != comm_size - 1 )
        {
            for ( int i = 0; i < 5; ++i )
            {
                parallel_search.getDomainEntitiesFromRange( point_ids[i],
                                                            domain_entities );
                TEST_EQUALITY( 1, domain_entities.size() );
                TEST_EQUALITY( point_ids[i], domain_entities[0] );
            }
        }
        if ( comm_rank != 0 )
        {
            for ( int i = 0; i < 5; ++i )
            {
                parallel_search.rangeParametricCoordinatesInDomain(
                    local_domain[i], local_range[i], range_coords );
                TEST_EQUALITY( range_coords[0], 0.5 );
                TEST_EQUALITY( range_coords[1], 0.5 );
                TEST_EQUALITY( range_coords[2], local_range[i] + 0.5 );
                TEST_EQUALITY(
                    comm_size - comm_rank - 1,
                    parallel_search.rangeEntityOwnerRank( local_range[i] ) );
            }
        }
    }

    // Check that proc zero had all points not found.
    int num_missed = ( comm_rank != comm_size - 1 ) ? 0 : 5;
    Teuchos::Array<EntityId> missed_ids(
        parallel_search.getMissedRangeEntityIds() );
    TEST_EQUALITY( missed_ids.size(), num_missed );
    std::sort( point_ids.begin(), point_ids.end() );
    std::sort( missed_ids.begin(), missed_ids.end() );
    for ( int i = 0; i < num_missed; ++i )
    {
        TEST_EQUALITY( missed_ids[i], point_ids[i] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, no_range_0_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes * ( comm_size - comm_rank - 1 ) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
            ->addEntity( BoxGeometry( id, comm_rank, id, 0.0, 0.0, id, 1.0, 1.0,
                                      id + 1.0 ) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_points = ( comm_rank != 0 ) ? 5 : 0;
    Teuchos::Array<double> point( 3 );
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points * comm_rank + i;
        point_ids[i] = id;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
            ->addEntity( Point( id, comm_rank, point ) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    Teuchos::Array<EntityId> local_range;
    Teuchos::Array<EntityId> local_domain;
    Teuchos::Array<EntityId> range_entities;
    if ( comm_rank != comm_size - 1 )
    {
        for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
              ++domain_it )
        {
            parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                        range_entities );
            TEST_EQUALITY( 1, range_entities.size() );
            TEST_EQUALITY( range_entities[0], domain_it->id() );
            local_range.push_back( range_entities[0] );
            local_domain.push_back( domain_it->id() );
        }
    }
    int range_size = ( comm_rank != comm_size - 1 ) ? 5 : 0;
    TEST_EQUALITY( local_range.size(), range_size );

    Teuchos::ArrayView<const double> range_coords;
    Teuchos::Array<EntityId> domain_entities;
    for ( int i = 0; i < num_points; ++i )
    {
        parallel_search.getDomainEntitiesFromRange( point_ids[i],
                                                    domain_entities );
        TEST_EQUALITY( 1, domain_entities.size() );
        TEST_EQUALITY( point_ids[i], domain_entities[0] );
    }

    for ( int i = 0; i < range_size; ++i )
    {
        parallel_search.rangeParametricCoordinatesInDomain(
            local_domain[i], local_range[i], range_coords );
        TEST_EQUALITY( range_coords[0], 0.5 );
        TEST_EQUALITY( range_coords[1], 0.5 );
        TEST_EQUALITY( range_coords[2], local_range[i] + 0.5 );
        TEST_EQUALITY( comm_size - comm_rank - 1,
                       parallel_search.rangeEntityOwnerRank( local_range[i] ) );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( parallel_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, many_to_many_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes * ( comm_size - comm_rank - 1 ) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
            ->addEntity( BoxGeometry( id, comm_rank, id, 0.0, 0.0, id, 1.0, 1.0,
                                      id + 1.0 ) );
    }

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_points = 10;
    Teuchos::Array<double> point( 3 );
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points * comm_rank + i;
        point_ids[i] = id;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = comm_rank * 5.0 + i + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
            ->addEntity( Point( id, comm_rank, point ) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    if ( comm_rank < comm_size - 1 )
    {
        Teuchos::Array<EntityId> local_range;
        Teuchos::Array<EntityId> local_domain;
        Teuchos::Array<EntityId> range_entities;
        Teuchos::Array<int> range_ranks;
        for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
              ++domain_it )
        {
            parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                        range_entities );
            TEST_EQUALITY( 2, range_entities.size() );
            local_range.push_back( range_entities[0] );
            local_domain.push_back( domain_it->id() );
            range_ranks.push_back(
                parallel_search.rangeEntityOwnerRank( range_entities[0] ) );
            local_range.push_back( range_entities[1] );
            local_domain.push_back( domain_it->id() );
            range_ranks.push_back(
                parallel_search.rangeEntityOwnerRank( range_entities[1] ) );
        }
        TEST_EQUALITY( local_range.size(), 10 );

        Teuchos::ArrayView<const double> range_coords;
        for ( int i = 0; i < 10; ++i )
        {
            parallel_search.rangeParametricCoordinatesInDomain(
                local_domain[i], local_range[i], range_coords );
            TEST_EQUALITY( range_coords[0], 0.5 );
            TEST_EQUALITY( range_coords[1], 0.5 );
            TEST_EQUALITY( range_coords[2],
                           local_range[i] - 5.0 * range_ranks[i] + 0.5 );
        }
    }
    else
    {
        Teuchos::Array<EntityId> local_range;
        Teuchos::Array<EntityId> local_domain;
        Teuchos::Array<EntityId> range_entities;
        Teuchos::Array<int> range_ranks;
        for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
              ++domain_it )
        {
            parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                        range_entities );
            TEST_EQUALITY( 1, range_entities.size() );
            TEST_EQUALITY( range_entities[0], domain_it->id() );
            local_range.push_back( range_entities[0] );
            local_domain.push_back( domain_it->id() );
            range_ranks.push_back(
                parallel_search.rangeEntityOwnerRank( range_entities[0] ) );
        }
        TEST_EQUALITY( local_range.size(), 5 );

        Teuchos::ArrayView<const double> range_coords;
        for ( int i = 0; i < 5; ++i )
        {
            parallel_search.rangeParametricCoordinatesInDomain(
                local_domain[i], local_range[i], range_coords );
            TEST_EQUALITY( range_coords[0], 0.5 );
            TEST_EQUALITY( range_coords[1], 0.5 );
            TEST_EQUALITY( range_coords[2],
                           local_range[i] - 5.0 * range_ranks.back() + 0.5 );
        }
        for ( int i = 0; i < 5; ++i )
        {
            TEST_EQUALITY( range_ranks[i], 0.0 );
        }
    }

    // Check that proc zero had some points not found.
    int num_missed = ( comm_rank != comm_size - 1 ) ? 0 : 5;
    Teuchos::Array<EntityId> missed_ids(
        parallel_search.getMissedRangeEntityIds() );
    TEST_EQUALITY( missed_ids.size(), num_missed );
    std::sort( point_ids.begin(), point_ids.end() );
    std::sort( missed_ids.begin(), missed_ids.end() );
    for ( int i = 0; i < num_missed; ++i )
    {
        TEST_EQUALITY( missed_ids[i], point_ids[i + 5] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, point_multiple_neighbors_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int box_id = comm_size - comm_rank - 1;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
        ->addEntity( BoxGeometry( box_id, box_id, box_id, 0.0, 0.0, box_id, 1.0,
                                  1.0, box_id + 1.0 ) );

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    Teuchos::Array<double> point( 3 );
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = comm_rank;
    int point_id = comm_rank;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
        ->addEntity( Point( point_id, comm_rank, point ) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    if ( comm_rank > 0 )
    {
        Teuchos::Array<EntityId> local_range;
        parallel_search.getRangeEntitiesFromDomain( box_id, local_range );
        TEST_EQUALITY( 2, local_range.size() );
        std::sort( local_range.begin(), local_range.end() );
        TEST_EQUALITY( Teuchos::as<int>( local_range[0] ),
                       comm_size - comm_rank - 1 );
        TEST_EQUALITY( Teuchos::as<int>( local_range[1] ),
                       comm_size - comm_rank );
        TEST_EQUALITY( parallel_search.rangeEntityOwnerRank( local_range[0] ),
                       comm_size - comm_rank - 1 );
        TEST_EQUALITY( parallel_search.rangeEntityOwnerRank( local_range[1] ),
                       comm_size - comm_rank );
        Teuchos::ArrayView<const double> range_centroid;
        parallel_search.rangeParametricCoordinatesInDomain(
            box_id, local_range[0], range_centroid );
        TEST_EQUALITY( range_centroid[0], 0.5 );
        TEST_EQUALITY( range_centroid[1], 0.5 );
        TEST_EQUALITY( range_centroid[2], comm_size - comm_rank - 1 );
        parallel_search.rangeParametricCoordinatesInDomain(
            box_id, local_range[1], range_centroid );
        TEST_EQUALITY( range_centroid[0], 0.5 );
        TEST_EQUALITY( range_centroid[1], 0.5 );
        TEST_EQUALITY( range_centroid[2], comm_size - comm_rank );
    }
    else
    {
        Teuchos::Array<EntityId> local_range;
        parallel_search.getRangeEntitiesFromDomain( box_id, local_range );
        TEST_EQUALITY( 1, local_range.size() );
        TEST_EQUALITY( Teuchos::as<int>( local_range[0] ),
                       comm_size - comm_rank - 1 );
        TEST_EQUALITY( parallel_search.rangeEntityOwnerRank( local_range[0] ),
                       comm_size - comm_rank - 1 );
        Teuchos::ArrayView<const double> range_centroid;
        parallel_search.rangeParametricCoordinatesInDomain(
            box_id, local_range[0], range_centroid );
        TEST_EQUALITY( range_centroid[0], 0.5 );
        TEST_EQUALITY( range_centroid[1], 0.5 );
        TEST_EQUALITY( range_centroid[2], comm_size - comm_rank - 1 );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( parallel_search.getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, global_missed_range_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes * ( comm_size - comm_rank - 1 ) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
            ->addEntity( BoxGeometry( id, comm_rank, id, 0.0, 0.0, id, 1.0, 1.0,
                                      id + 1.0 ) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_points = 5;
    Teuchos::Array<double> point( 3 );
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points + 1 );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points * comm_rank + i;
        point_ids[i] = id;
        point[0] = 0.5;
        point[1] = 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
            ->addEntity( Point( id, comm_rank, point ) );
    }

    // Add a bad point.
    id = num_points * comm_rank + 1000;
    point_ids[5] = id;
    point[0] = -100.0;
    point[1] = 0.0;
    point[2] = 0.0;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
        ->addEntity( Point( id, comm_rank, point ) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    Teuchos::Array<EntityId> local_range;
    Teuchos::Array<EntityId> local_domain;
    Teuchos::Array<EntityId> range_entities;
    for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
          ++domain_it )
    {
        parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                    range_entities );
        TEST_EQUALITY( 1, range_entities.size() );
        TEST_EQUALITY( range_entities[0], domain_it->id() );
        local_range.push_back( range_entities[0] );
        local_domain.push_back( domain_it->id() );
    }
    TEST_EQUALITY( local_range.size(), 5 );
    Teuchos::ArrayView<const double> range_coords;
    Teuchos::Array<EntityId> domain_entities;
    for ( int i = 0; i < 5; ++i )
    {
        parallel_search.getDomainEntitiesFromRange( point_ids[i],
                                                    domain_entities );
        TEST_EQUALITY( 1, domain_entities.size() );
        TEST_EQUALITY( point_ids[i], domain_entities[0] );
    }
    for ( int i = 0; i < 5; ++i )
    {
        parallel_search.rangeParametricCoordinatesInDomain(
            local_domain[i], local_range[i], range_coords );
        TEST_EQUALITY( range_coords[0], 0.5 );
        TEST_EQUALITY( range_coords[1], 0.5 );
        TEST_EQUALITY( range_coords[2], local_range[i] + 0.5 );
        TEST_EQUALITY( comm_size - comm_rank - 1,
                       parallel_search.rangeEntityOwnerRank( local_range[i] ) );
    }

    // Check that the bad point was found.
    Teuchos::ArrayView<const EntityId> missed_range =
        parallel_search.getMissedRangeEntityIds();
    TEST_EQUALITY( missed_range.size(), 1 );
    TEST_EQUALITY( missed_range[0],
                   Teuchos::as<EntityId>( num_points * comm_rank + 1000 ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, local_missed_range_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_boxes = 5;
    int id = 0;
    for ( int i = 0; i < num_boxes; ++i )
    {
        id = num_boxes * ( comm_size - comm_rank - 1 ) + i;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( domain_set )
            ->addEntity( BoxGeometry( id, comm_rank, id, id, id, id, id + 1.0,
                                      id + 1.0, id + 1.0 ) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( 3 );

    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    plist.set<bool>( "Track Missed Range Entities", true );
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
        Teuchos::rcp( new BasicEntitySet( comm, 3 ) );
    int num_points = 5;
    Teuchos::Array<double> point( 3 );
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points + 1 );
    for ( int i = 0; i < num_points; ++i )
    {
        id = num_points * comm_rank + i;
        point_ids[i] = id;
        point[0] = id + 0.5;
        point[1] = id + 0.5;
        point[2] = id + 0.5;
        Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
            ->addEntity( Point( id, comm_rank, point ) );
    }

    // Add a bad point.
    id = num_points * comm_rank;
    point_ids[5] = id + 1000;
    point[0] = id + 0.5;
    point[1] = id + 0.5;
    point[2] = id + 1.5;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>( range_set )
        ->addEntity( Point( point_ids[5], comm_rank, point ) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map =
        Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( 0 );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    Teuchos::Array<EntityId> local_range;
    Teuchos::Array<EntityId> local_domain;
    Teuchos::Array<EntityId> range_entities;
    for ( domain_it = domain_it.begin(); domain_it != domain_it.end();
          ++domain_it )
    {
        parallel_search.getRangeEntitiesFromDomain( domain_it->id(),
                                                    range_entities );
        TEST_EQUALITY( 1, range_entities.size() );
        TEST_EQUALITY( range_entities[0], domain_it->id() );
        local_range.push_back( range_entities[0] );
        local_domain.push_back( domain_it->id() );
    }
    TEST_EQUALITY( local_range.size(), 5 );
    Teuchos::ArrayView<const double> range_coords;
    Teuchos::Array<EntityId> domain_entities;
    for ( int i = 0; i < 5; ++i )
    {
        parallel_search.getDomainEntitiesFromRange( point_ids[i],
                                                    domain_entities );
        TEST_EQUALITY( 1, domain_entities.size() );
        TEST_EQUALITY( point_ids[i], domain_entities[0] );
    }
    for ( int i = 0; i < 5; ++i )
    {
        parallel_search.rangeParametricCoordinatesInDomain(
            local_domain[i], local_range[i], range_coords );
        TEST_EQUALITY( range_coords[0], local_range[i] + 0.5 );
        TEST_EQUALITY( range_coords[1], local_range[i] + 0.5 );
        TEST_EQUALITY( range_coords[2], local_range[i] + 0.5 );
        TEST_EQUALITY( comm_size - comm_rank - 1,
                       parallel_search.rangeEntityOwnerRank( local_range[i] ) );
    }

    // Check that the bad point was found.
    Teuchos::ArrayView<const EntityId> missed_range =
        parallel_search.getMissedRangeEntityIds();
    TEST_EQUALITY( missed_range.size(), 1 );
    TEST_EQUALITY( missed_range[0],
                   Teuchos::as<EntityId>( num_points * comm_rank + 1000 ) );
}

//---------------------------------------------------------------------------//
// end tstParallelSearch.cpp
//---------------------------------------------------------------------------//
