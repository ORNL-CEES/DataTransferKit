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
#include <DTK_Entity.hpp>
#include <DTK_Box.hpp>
#include <DTK_AbstractObjectRegistry.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_AbstractFactoryStd.hpp>

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
// Local point set test.
TEUCHOS_UNIT_TEST( Point, set_test )
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
    Entity point_1 = Point<3>(0, comm_rank, p1);

    // Make a second point.
    double x_2 = 3.2 - comm_rank;
    double y_2 = -9.233 - comm_rank;
    double z_2 = 1.3 - comm_rank;
    Teuchos::Array<double> p2(3);
    p2[0] = x_2;
    p2[1] = y_2;
    p2[2] = z_2;
    Entity point_2 = Point<3>(1, comm_rank, p2);

    // Make an entity set.
    Teuchos::RCP<EntitySet> entity_set = Teuchos::rcp(
	new BasicEntitySet(comm, 3) );

    // Add the points to the set.
    entity_set->addEntity( point_1 );
    entity_set->addEntity( point_2 );

    // Get an iterator to the entity set objects.
    AbstractIterator<Entity> node_it =
	entity_set->entityIterator( NODE );
    AbstractIterator<Entity> edge_it =
	entity_set->entityIterator( EDGE );
    AbstractIterator<Entity> face_it =
	entity_set->entityIterator( FACE );
    AbstractIterator<Entity> volume_it =
	entity_set->entityIterator( VOLUME );

    // Check the entity set.
    TEST_EQUALITY( entity_set->name(), "DTK Basic Entity Set" );
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );
    TEST_EQUALITY( node_it.size(), 2 );
    TEST_EQUALITY( edge_it.size(), 0 );
    TEST_EQUALITY( face_it.size(), 0 );
    TEST_EQUALITY( volume_it.size(), 0 );

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
    entity_set->getEntity( 0, entity );
    TEST_EQUALITY( 0, entity.id() );
    entity_set->getEntity( 1, entity );
    TEST_EQUALITY( 1, entity.id() );

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
    TEST_FLOATING_EQUALITY( global_bounds[0], 3.2 - comm_size + 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[1], -9.233 - comm_size + 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[2], 1.3 - comm_size + 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[3], 3.2 + comm_size - 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[4], -9.233 + comm_size - 1.0, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[5], 1.3 + comm_size - 1.0, 1.0e-12 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Point, modification_test )
{
    using namespace DataTransferKit;

    // Register the Entity classes.
    AbstractObjectRegistry<Entity,Point<3> >::registerDerivedClasses();

    // Register the EntitySet classes.
    AbstractObjectRegistry<EntitySet,BasicEntitySet>::registerDerivedClasses();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Create a builder for the entity sets.
    Teuchos::RCP<AbstractBuilder<EntitySet> > builder = EntitySet::getBuilder();

    // Make an entity set on process 0.
    Teuchos::RCP<EntitySet> entity_set;
    Teuchos::Array<Entity> points(2);
    int entity_set_key = -1;
    double x_1 = 3.2 + comm_size;
    double y_1 = -9.233 + comm_size;
    double z_1 = 1.3 + comm_size;
    double x_2 = 3.2 - comm_size;
    double y_2 = -9.233 - comm_size;
    double z_2 = 1.3 - comm_size;
    if ( 0 == comm->getRank() )
    {
	Teuchos::Array<double> p1(3);
	p1[0] = x_1;
	p1[1] = y_1;
	p1[2] = z_1;
	points[0] = Point<3>(0, 0, p1);
	Teuchos::Array<double> p2(3);
	p2[0] = x_2;
	p2[1] = y_2;
	p2[2] = z_2;
	points[1] = Point<3>(1, 0, p2);

	entity_set = Teuchos::rcp(new BasicEntitySet(comm,3) );
	entity_set_key = builder->getIntegralKey( entity_set->name() );
    }

    // Create an entity set.
    Teuchos::broadcast( *comm, 0, Teuchos::Ptr<int>(&entity_set_key) );
    entity_set = builder->create( entity_set_key );
    entity_set->assignCommunicator( comm );
    TEST_EQUALITY( entity_set->physicalDimension(), 0 );

    // Broadcast the points with indirect serialization through the geometric
    // entity api.
    Teuchos::broadcast( *comm, 0, points() );

    // Add the points to the entity set.
    entity_set->addEntity( points[0] );
    entity_set->addEntity( points[1] );
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );

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
    TEST_FLOATING_EQUALITY( global_bounds[0], x_2, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[1], y_2, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[2], z_2, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[3], x_1, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[4], y_1, 1.0e-12 );
    TEST_FLOATING_EQUALITY( global_bounds[5], z_1, 1.0e-12 );
}

//---------------------------------------------------------------------------//
// end tstPoint.cpp
//---------------------------------------------------------------------------//

