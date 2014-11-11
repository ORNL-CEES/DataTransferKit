//---------------------------------------------------------------------------//
/*!
 * \file tstBasicEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief BasicEntitySet unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_BasicEntitySet.hpp>
#include <DTK_Point.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( BasicEntitySet, basic_entity_set_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make point.
    double x_1 = 3.2 + comm_rank;
    double y_1 = -9.233 + comm_rank;
    double z_1 = 1.3 + comm_rank;
    Teuchos::Array<double> p1(3);
    p1[0] = x_1;
    p1[1] = y_1;
    p1[2] = z_1;
    Entity point_1 = Point(0, comm_rank, p1);

    // Make a second point.
    double x_2 = 3.2 - comm_rank;
    double y_2 = -9.233 - comm_rank;
    double z_2 = 1.3 - comm_rank;
    Teuchos::Array<double> p2(3);
    p2[0] = x_2;
    p2[1] = y_2;
    p2[2] = z_2;
    Entity point_2 = Point(1, comm_rank, p2);

    // Make a box.
    Entity box_1 = Box( 2, comm_rank, 1,
			std::min(p1[0],p2[0]),
			std::min(p1[1],p2[1]),
			std::min(p1[2],p2[2]),
			std::max(p1[0],p2[0]),
			std::max(p1[1],p2[1]),
			std::max(p1[2],p2[2]) );

    // Make an entity set.
    Teuchos::RCP<EntitySet> entity_set = 
	Teuchos::rcp( new BasicEntitySet(comm, 3) );

    // Add the points and box to the set.
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(entity_set)->addEntity( point_1 );
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(entity_set)->addEntity( point_2 );
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(entity_set)->addEntity( box_1 );

    // Get an iterator to the entity set objects.
    EntityIterator node_it = entity_set->entityIterator( ENTITY_TYPE_NODE );
    EntityIterator edge_it = entity_set->entityIterator( ENTITY_TYPE_EDGE );
    EntityIterator face_it = entity_set->entityIterator( ENTITY_TYPE_FACE );
    EntityIterator volume_it = entity_set->entityIterator( ENTITY_TYPE_VOLUME );

    // Check the entity set.
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );
    TEST_EQUALITY( node_it.size(), 2 );
    TEST_EQUALITY( edge_it.size(), 0 );
    TEST_EQUALITY( face_it.size(), 0 );
    TEST_EQUALITY( volume_it.size(), 1 );

    // Check the nodes.
    node_it = node_it.begin();
    Entity ge0 = *node_it;
    ++node_it;
    Entity ge1 = *node_it;
    if ( ge1.id() == 1 )
    {
	TEST_EQUALITY( ge0.id(), 0 );
    }
    else if ( ge1.id() == 0 )
    {
	TEST_EQUALITY( ge0.id(), 1 );
    }
    else
    {
	TEST_ASSERT( false );
    }
    Entity entity;
    entity_set->getEntity( ENTITY_TYPE_NODE, 0, entity );
    TEST_EQUALITY( 0, entity.id() );
    entity_set->getEntity( ENTITY_TYPE_NODE, 1, entity );
    TEST_EQUALITY( 1, entity.id() );

    // Check the box.
    Entity ge2 = *volume_it;
    TEST_EQUALITY( ge2.id(), 2 );
    TEST_EQUALITY( ge2.ownerRank(), comm_rank );
    TEST_ASSERT( !ge2.inBlock(0) );
    TEST_ASSERT( ge2.inBlock(1) );

    // Check the bounding boxes.
    Teuchos::Tuple<double,6> local_bounds;
    entity_set->localBoundingBox( local_bounds );
    TEST_EQUALITY( local_bounds[0], x_2 );
    TEST_EQUALITY( local_bounds[1], y_2 );
    TEST_EQUALITY( local_bounds[2], z_2 );
    TEST_EQUALITY( local_bounds[3], x_1 );
    TEST_EQUALITY( local_bounds[4], y_1 );
    TEST_EQUALITY( local_bounds[5], z_1 );

    Teuchos::Tuple<double,6> global_bounds;
    entity_set->globalBoundingBox( global_bounds );
    TEST_FLOATING_EQUALITY( global_bounds[0], 3.2 - (comm_size-1.0), 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[1], -9.233 - (comm_size-1.0), 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[2], 1.3 - (comm_size-1.0), 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[3], 3.2 + (comm_size-1.0), 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[4], -9.233 + (comm_size-1.0), 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[5], 1.3 + (comm_size-1.0), 1.0e-12 );
}

//---------------------------------------------------------------------------//
// end tstBasicEntitySet.cpp
//---------------------------------------------------------------------------//

