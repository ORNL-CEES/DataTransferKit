//---------------------------------------------------------------------------//
/*! 
 * \file tstParallelSearch.cpp
 * \author Stuart R. Slattery
 * \brief ParallelSearch unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_ParallelSearch.hpp>
#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_Box.hpp>
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
TEUCHOS_UNIT_TEST( ParallelSearch, all_to_all_test )
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
    for ( int i = 0; i < num_boxes; ++i )
    {
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(i,comm_rank,i,0.0,0.0,i,1.0,1.0,i+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( ENTITY_TYPE_VOLUME );
    
    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

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
	bool on_surface = (i%2==0);
	int id = num_points*comm_rank + i;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(id,comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( ENTITY_TYPE_NODE );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    Teuchos::Array<EntityId> local_range;
    Teuchos::Array<EntityId> range_entities;
    for ( domain_it = domain_it.begin();
	  domain_it != domain_it.end();
	  ++domain_it )
    {
	parallel_search.getRangeEntitiesFromDomain(
	    domain_it->id(), range_entities );
	TEST_EQUALITY( comm_size, range_entities.size() );
	for ( int n = 0; n < comm_size; ++n )
	{
	    TEST_EQUALITY( range_entities[n] % num_points, domain_it->id() );
	    local_range.push_back( range_entities[n] );
	}
    }
    TEST_EQUALITY( local_range.size(), num_points*comm_size );
    Teuchos::ArrayView<const double> range_coords;
    Teuchos::Array<EntityId> domain_entities;
    for ( int i = 0; i < num_points*comm_size; ++i )
    {
	parallel_search.getDomainEntitiesFromRange(
	    local_range[i], domain_entities );
	TEST_EQUALITY( 1, domain_entities.size() );
	TEST_EQUALITY( local_range[i] % num_points, domain_entities[0] );
	parallel_search.rangeParametricCoordinatesInDomain( 
	    domain_entities[0], local_range[i], range_coords );
	TEST_EQUALITY( range_coords[0], 0.5 );
	TEST_EQUALITY( range_coords[1], 0.5 );
	TEST_EQUALITY( range_coords[2], domain_entities[0] + 0.5 );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, one_to_one_test )
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
	    Box(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( ENTITY_TYPE_VOLUME );
    
    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    for ( int i = 0; i < num_points; ++i )
    {
	bool on_surface = (i%2==0);
	id = num_points*comm_rank + i;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = id + 0.5;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(id,comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( ENTITY_TYPE_NODE );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    Teuchos::Array<EntityId> local_range;
    Teuchos::Array<EntityId> range_entities;
    for ( domain_it = domain_it.begin();
	  domain_it != domain_it.end();
	  ++domain_it )
    {
	parallel_search.getRangeEntitiesFromDomain(
	    domain_it->id(), range_entities );
	TEST_EQUALITY( 1, range_entities.size() );
	TEST_EQUALITY( range_entities[0], domain_it->id() );
	local_range.push_back( range_entities[0] );
    }
    TEST_EQUALITY( local_range.size(), 5 );
    Teuchos::ArrayView<const double> range_coords;
    Teuchos::Array<EntityId> domain_entities;
    for ( int i = 0; i < 5; ++i )
    {
	parallel_search.getDomainEntitiesFromRange(
	    local_range[i], domain_entities );
	TEST_EQUALITY( 1, domain_entities.size() );
	TEST_EQUALITY( local_range[i], domain_entities[0] );
	parallel_search.rangeParametricCoordinatesInDomain( 
	    domain_entities[0], local_range[i], range_coords );
	TEST_EQUALITY( range_coords[0], 0.5 );
	TEST_EQUALITY( range_coords[1], 0.5 );
	TEST_EQUALITY( range_coords[2], local_range[i] + 0.5 );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ParallelSearch, no_domain_0_test )
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
		Box(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
	}
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over all of the boxes.
    EntityIterator domain_it = domain_set->entityIterator( ENTITY_TYPE_VOLUME );
    
    // Build a parallel search over the boxes.
    Teuchos::ParameterList plist;
    ParallelSearch parallel_search( comm, 3, domain_it, domain_map, plist );

    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    for ( int i = 0; i < num_points; ++i )
    {
	bool on_surface = (i%2==0);
	id = num_points*comm_rank + i;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = id + 0.5;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(id,comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Get an iterator over the points.
    EntityIterator range_it = range_set->entityIterator( ENTITY_TYPE_NODE );

    // Do the search.
    parallel_search.search( range_it, range_map, plist );

    // Check the results of the search.
    if ( comm_rank > 0 )
    {
	Teuchos::Array<EntityId> local_range;
	Teuchos::Array<EntityId> range_entities;
	for ( domain_it = domain_it.begin();
	      domain_it != domain_it.end();
	      ++domain_it )
	{
	    parallel_search.getRangeEntitiesFromDomain(
		domain_it->id(), range_entities );
	    TEST_EQUALITY( 1, range_entities.size() );
	    TEST_EQUALITY( range_entities[0], domain_it->id() );
	    local_range.push_back( range_entities[0] );
	}
	TEST_EQUALITY( local_range.size(), 5 );
	Teuchos::ArrayView<const double> range_coords;
	Teuchos::Array<EntityId> domain_entities;
	for ( int i = 0; i < 5; ++i )
	{
	    parallel_search.getDomainEntitiesFromRange(
		local_range[i], domain_entities );
	    TEST_EQUALITY( 1, domain_entities.size() );
	    TEST_EQUALITY( local_range[i], domain_entities[0] );
	    parallel_search.rangeParametricCoordinatesInDomain( 
		domain_entities[0], local_range[i], range_coords );
	    TEST_EQUALITY( range_coords[0], 0.5 );
	    TEST_EQUALITY( range_coords[1], 0.5 );
	    TEST_EQUALITY( range_coords[2], local_range[i] + 0.5 );
	}
    }
}

// //---------------------------------------------------------------------------//
// TEUCHOS_UNIT_TEST( ParallelSearch, no_range_0_test )
// {
//     using namespace DataTransferKit;

//     // Get the communicator.
//     Teuchos::RCP<const Teuchos::Comm<int> > comm =
// 	Teuchos::DefaultComm<int>::getComm();
//     int comm_rank = comm->getRank();
//     int comm_size = comm->getSize();

//     // Make a domain entity set.
//     Teuchos::RCP<EntitySet> domain_set = 
// 	Teuchos::rcp( new BasicEntitySet(comm,3) );
//     int num_boxes = 5;
//     int id = 0;
//     for ( int i = 0; i < num_boxes; ++i )
//     {
// 	id = num_boxes*(comm_size-comm_rank-1) + i;
// 	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
// 	    Box(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
//     }

//     // Get an iterator over all of the boxes.
//     EntityIterator domain_it = domain_set->entityIterator( ENTITY_TYPE_VOLUME );
    
//     // Build a parallel search over the boxes.
//     Teuchos::ParameterList plist;
//     ParallelSearch parallel_search( comm, 3, domain_it, plist );

//     // Make a range entity set.
//     Teuchos::RCP<EntitySet> range_set =
// 	Teuchos::rcp( new BasicEntitySet(comm,3) );
//     int num_points = 5;
//     if ( comm_rank > 0 )
//     {
// 	Teuchos::Array<double> point(3);
// 	for ( int i = 0; i < num_points; ++i )
// 	{
// 	    bool on_surface = (i%2==0);
// 	    id = num_points*comm_rank + i;
// 	    point[0] = 0.5;
// 	    point[1] = 0.5;
// 	    point[2] = id + 0.5;
// 	    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
// 		Point(id,comm_rank,point,on_surface) );
// 	}
//     }

//     // Construct a local map for the points.
//     Teuchos::RCP<EntityLocalMap> range_map = 
// 	Teuchos::rcp( new BasicGeometryLocalMap() );

//     // Get an iterator over the points.
//     EntityIterator range_it = range_set->entityIterator( ENTITY_TYPE_NODE );

//     // Search the tree.
//     Teuchos::Array<EntityId> range_ids;
//     Teuchos::Array<int> range_ranks;
//     Teuchos::Array<double> range_centroids;
//     parallel_search.search( range_it, range_map, plist,
// 				 range_ids, range_ranks, range_centroids );

//     // Check the results of the search.
//     if ( comm_rank < comm_size - 1 )
//     {
// 	int num_recv = num_points;
// 	TEST_EQUALITY( num_recv, range_ids.size() );
// 	TEST_EQUALITY( num_recv, range_ranks.size() );
// 	TEST_EQUALITY( 3*num_recv, range_centroids.size() );
// 	for ( int i = 0; i < num_recv; ++i )
// 	{
// 	    TEST_EQUALITY( range_ranks[i], comm_size - comm_rank - 1 );
// 	    TEST_EQUALITY( range_centroids[3*i], 0.5 );
// 	    TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
// 	    TEST_EQUALITY( range_centroids[3*i+2], range_ids[i]+0.5 );
// 	}
// 	std::sort( range_ids.begin(), range_ids.end() );
// 	for ( int j = 0; j < num_recv; ++j )
// 	{
// 	    id = num_points*(comm_size-comm_rank-1) + j;
// 	    TEST_EQUALITY( Teuchos::as<int>(range_ids[j]), id );
// 	}
//     }
//     else
//     {
// 	TEST_EQUALITY( 0, range_ids.size() );
// 	TEST_EQUALITY( 0, range_ranks.size() );
// 	TEST_EQUALITY( 0, range_centroids.size() );
//     }
// }

// //---------------------------------------------------------------------------//
// TEUCHOS_UNIT_TEST( ParallelSearch, many_to_many_test )
// {
//     using namespace DataTransferKit;

//     // Get the communicator.
//     Teuchos::RCP<const Teuchos::Comm<int> > comm =
// 	Teuchos::DefaultComm<int>::getComm();
//     int comm_rank = comm->getRank();
//     int comm_size = comm->getSize();

//     // Make a domain entity set.
//     Teuchos::RCP<EntitySet> domain_set = 
// 	Teuchos::rcp( new BasicEntitySet(comm,3) );
//     int num_boxes = 5;
//     int id = 0;
//     for ( int i = 0; i < num_boxes; ++i )
//     {
// 	id = num_boxes*(comm_size-comm_rank-1) + i;
// 	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
// 	    Box(id,comm_rank,id,0.0,0.0,id,1.0,1.0,id+1.0) );
//     }

//     // Get an iterator over all of the boxes.
//     EntityIterator domain_it = domain_set->entityIterator( ENTITY_TYPE_VOLUME );
    
//     // Build a parallel search over the boxes.
//     Teuchos::ParameterList plist;
//     ParallelSearch parallel_search( comm, 3, domain_it, plist );

//     // Make a range entity set.
//     Teuchos::RCP<EntitySet> range_set =
// 	Teuchos::rcp( new BasicEntitySet(comm,3) );
//     int num_points = 10;
//     Teuchos::Array<double> point(3);
//     for ( int i = 0; i < num_points; ++i )
//     {
// 	bool on_surface = (i%2==0);
// 	id = num_points*comm_rank + i;
// 	point[0] = 0.5;
// 	point[1] = 0.5;
// 	point[2] = comm_rank*5.0 + i + 0.5;
// 	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
// 	    Point(id,comm_rank,point,on_surface) );
//     }

//     // Construct a local map for the points.
//     Teuchos::RCP<EntityLocalMap> range_map = 
// 	Teuchos::rcp( new BasicGeometryLocalMap() );

//     // Get an iterator over the points.
//     EntityIterator range_it = range_set->entityIterator( ENTITY_TYPE_NODE );

//     // Search the tree.
//     Teuchos::Array<EntityId> range_ids;
//     Teuchos::Array<int> range_ranks;
//     Teuchos::Array<double> range_centroids;
//     parallel_search.search( range_it, range_map, plist,
// 				 range_ids, range_ranks, range_centroids );

//     // Check the results of the search.
//     if ( comm_rank < comm_size - 1 )
//     {
// 	int num_recv = num_points;
// 	TEST_EQUALITY( num_recv, range_ids.size() );
// 	TEST_EQUALITY( num_recv, range_ranks.size() );
// 	TEST_EQUALITY( 3*num_recv, range_centroids.size() );
// 	for ( int i = 0; i < num_recv; ++i )
// 	{
// 	    TEST_EQUALITY( range_centroids[3*i], 0.5 );
// 	    TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
// 	    TEST_EQUALITY( range_centroids[3*i+2], 
// 			   range_ids[i]-5.0*range_ranks[i]+0.5 );
// 	}
// 	std::sort( range_ranks.begin(), range_ranks.end() );
// 	int k = 0;
// 	for ( int j = 1; j > -1; --j )
// 	{
// 	    for ( int i = 0; i < num_points / 2; ++i, ++k )
// 	    {
// 		TEST_EQUALITY( range_ranks[k], comm_size - comm_rank - 1 - j );
// 	    }
// 	}
//     }
//     else
//     {
// 	int num_recv = num_points / 2;
// 	TEST_EQUALITY( num_recv, range_ids.size() );
// 	TEST_EQUALITY( num_recv, range_ranks.size() );
// 	TEST_EQUALITY( 3*num_recv, range_centroids.size() );
// 	for ( int i = 0; i < num_recv; ++i )
// 	{
// 	    TEST_EQUALITY( range_centroids[3*i], 0.5 );
// 	    TEST_EQUALITY( range_centroids[3*i+1], 0.5 );
// 	    TEST_EQUALITY( range_centroids[3*i+2], 
// 			   range_ids[i]-5.0*range_ranks[i]+0.5 );
// 	    TEST_EQUALITY( range_ranks[i], 0 );
// 	}
//     }
// }

// //---------------------------------------------------------------------------//
// TEUCHOS_UNIT_TEST( ParallelSearch, point_multiple_neighbors_test )
// {
//     using namespace DataTransferKit;

//     // Get the communicator.
//     Teuchos::RCP<const Teuchos::Comm<int> > comm =
// 	Teuchos::DefaultComm<int>::getComm();
//     int comm_rank = comm->getRank();
//     int comm_size = comm->getSize();

//     // Make a domain entity set.
//     Teuchos::RCP<EntitySet> domain_set = 
// 	Teuchos::rcp( new BasicEntitySet(comm,3) );
//     int id = comm_size - comm_rank - 1;
//     Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
// 	Box(id,id,id,0.0,0.0,id,1.0,1.0,id+1.0) );

//     // Get an iterator over all of the boxes.
//     EntityIterator domain_it = domain_set->entityIterator( ENTITY_TYPE_VOLUME );
    
//     // Build a parallel search over the boxes.
//     Teuchos::ParameterList plist;
//     ParallelSearch parallel_search( comm, 3, domain_it, plist );

//     // Make a range entity set.
//     Teuchos::RCP<EntitySet> range_set =
// 	Teuchos::rcp( new BasicEntitySet(comm,3) );
//     Teuchos::Array<double> point(3);
//     point[0] = 0.5;
//     point[1] = 0.5;
//     point[2] = comm_rank;
//     id = comm_rank;
//     Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
// 	Point(id,comm_rank,point,false) );

//     // Construct a local map for the points.
//     Teuchos::RCP<EntityLocalMap> range_map = 
// 	Teuchos::rcp( new BasicGeometryLocalMap() );

//     // Get an iterator over the points.
//     EntityIterator range_it = range_set->entityIterator( ENTITY_TYPE_NODE );

//     // Search the tree.
//     Teuchos::Array<EntityId> range_ids;
//     Teuchos::Array<int> range_ranks;
//     Teuchos::Array<double> range_centroids;
//     parallel_search.search( range_it, range_map, plist,
//     				 range_ids, range_ranks, range_centroids );

//     // Check the results of the search.
//     if ( comm_rank > 0 )
//     {
//     	TEST_EQUALITY( 2, range_ids.size() );
//     	TEST_EQUALITY( 2, range_ranks.size() );
//     	TEST_EQUALITY( 6, range_centroids.size() );
//     	std::sort( range_ids.begin(), range_ids.end() );
//     	std::sort( range_ranks.begin(), range_ranks.end() );
//     	TEST_EQUALITY( Teuchos::as<int>(range_ids[0]), comm_size-comm_rank-2 );
//     	TEST_EQUALITY( Teuchos::as<int>(range_ids[1]), comm_size-comm_rank-1);
//     	TEST_EQUALITY( range_ranks[0], comm_size-comm_rank-2 );
//     	TEST_EQUALITY( range_ranks[1], comm_size-comm_rank-1 );
//     	TEST_EQUALITY( range_centroids[0], 0.5 );
//     	TEST_EQUALITY( range_centroids[1], 0.5 );
//     	TEST_EQUALITY( range_centroids[2], comm_size-comm_rank-2 );
//     	TEST_EQUALITY( range_centroids[3], 0.5 );
//     	TEST_EQUALITY( range_centroids[4], 0.5 );
//     	TEST_EQUALITY( range_centroids[5], comm_size-comm_rank-1 );
//     }
//     else
//     {
//     	TEST_EQUALITY( 1, range_ids.size() );
//     	TEST_EQUALITY( 1, range_ranks.size() );
//     	TEST_EQUALITY( 3, range_centroids.size() );
//     	TEST_EQUALITY( Teuchos::as<int>(range_ids[0]), comm_size-comm_rank-1 );
//     	TEST_EQUALITY( Teuchos::as<int>(range_ranks[0]), comm_size-comm_rank-1 );
//     	TEST_EQUALITY( range_centroids[0], 0.5 );
//     	TEST_EQUALITY( range_centroids[1], 0.5 );
//     	TEST_EQUALITY( range_centroids[2], comm_size-comm_rank-1 );
//     }
// }

//---------------------------------------------------------------------------//
// end tstParallelSearch.cpp
//---------------------------------------------------------------------------//
