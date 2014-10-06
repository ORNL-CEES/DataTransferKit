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
#include <DTK_GeometricEntity.hpp>
#include <DTK_Box.hpp>
#include <DTK_DerivedObjectRegistry.hpp>

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
    GeometricEntity point_1 = Point<3>(0, comm_rank, p1);

    // Make a second point.
    double x_2 = 3.2 - comm_rank;
    double y_2 = -9.233 - comm_rank;
    double z_2 = 1.3 - comm_rank;
    Teuchos::Array<double> p2(3);
    p2[0] = x_2;
    p2[1] = y_2;
    p2[2] = z_2;
    GeometricEntity point_2 = Point<3>(1, comm_rank, p2);

    // Make an entity set.
    Teuchos::RCP<EntitySet> entity_set = Teuchos::rcp(
	new BasicEntitySet(comm, 3) );

    // Add the points to the set.
    entity_set->addEntity( point_1 );
    entity_set->addEntity( point_2 );

    // Check the entity set.
    TEST_EQUALITY( entity_set->entitySetType(), "DTK Basic Entity Set" );
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );
    TEST_EQUALITY( entity_set->localNumberOfEntities(0), 2 );
    TEST_EQUALITY( entity_set->localNumberOfEntities(1), 0 );
    TEST_EQUALITY( entity_set->globalNumberOfEntities(0), 
		   Teuchos::as<std::size_t>(2*comm->getSize()) );
    TEST_EQUALITY( entity_set->globalNumberOfEntities(1), 0 );

    // Check the entities.
    Teuchos::Array<EntityId> ids( 2 );
    entity_set->localEntityIds( 0, ids );
    if ( ids[1] == 1 )
    {
	TEST_EQUALITY( ids[0], 0 );
    }
    else if ( ids[1] == 0 )
    {
	TEST_EQUALITY( ids[0], 1 );
    }
    else
    {
	TEST_ASSERT( false );
    }
    GeometricEntity entity;
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

    // Register the GeometricEntity classes.
    DerivedObjectRegistry<GeometricEntity,Point<3> >::registerDerivedClasses();

    // Register the EntitySet classes.
    DerivedObjectRegistry<EntitySet,BasicEntitySet>::registerDerivedClasses();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Create a builder for the entity sets.
    Teuchos::RCP<AbstractBuilder<EntitySet> > builder = EntitySet::getBuilder();

    // Make an entity set on process 0.
    Teuchos::RCP<EntitySet> entity_set;
    Teuchos::Array<GeometricEntity> points(2);
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
	entity_set_key = builder->getIntegralKey( entity_set->entitySetType() );
    }

    // Create an entity set.
    Teuchos::broadcast( *comm, 0, Teuchos::Ptr<int>(&entity_set_key) );
    entity_set = builder->create( entity_set_key );
    entity_set->assignCommunicator( comm );

    // Broadcast the points with indirect serialization through the geometric
    // entity api.
    Teuchos::broadcast( *comm, 0, points() );

    // Add the points to the entity set.
    entity_set->addEntity( points[0] );
    entity_set->addEntity( points[1] );

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

