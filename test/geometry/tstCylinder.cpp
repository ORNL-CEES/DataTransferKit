//---------------------------------------------------------------------------//
/*!
 * \file tstCylinder.cpp
 * \author Stuart R. Slattery
 * \brief Bounding Box unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_BoundingBox.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_GeometricEntity.hpp>

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
// Helper Functions
//---------------------------------------------------------------------------//
bool softEquivalence( double a1, double a2, double tol=1.0e-6 )
{
    if ( std::abs( a1 - a2 ) < tol )
    {
	return true;
    }
    else
    {
	return false;
    }
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Cylinder, cylinder_test )
{
    using namespace DataTransferKit;

    // Make sure that PI is PI.
    double zero = 0.0;
    double pi = 2.0 * std::acos(zero);
    TEST_ASSERT( softEquivalence( pi, 3.14159, 1.0e-5 ) );

    // Build a series of random cylinders.
    int num_cylinders = 100;
    for ( int i = 0; i < num_cylinders; ++i )
    {
	// Make a cylinder.
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5;
	Cylinder cylinder( length, radius, centroid_x, centroid_y, centroid_z );

	// Check the geometric bounds.
	TEST_ASSERT( cylinder.length() == length );
	TEST_ASSERT( cylinder.radius() == radius );
	TEST_ASSERT( cylinder.centroid()[0] == centroid_x );
	TEST_ASSERT( cylinder.centroid()[1] == centroid_y );
	TEST_ASSERT( cylinder.centroid()[2] == centroid_z );

	// Compute the measure.
	double measure = pi*radius*radius*length;
	TEST_ASSERT( softEquivalence( cylinder.measure(), measure, 1.0e-6 ) );

	// Compute the bounding box.
	BoundingBox box = cylinder.boundingBox();
	Teuchos::Tuple<double,6> box_bounds = box.getBounds();
	TEST_ASSERT( box_bounds[0] == centroid_x - radius );
	TEST_ASSERT( box_bounds[1] == centroid_y - radius );
	TEST_ASSERT( box_bounds[2] == centroid_z - length/2 );
	TEST_ASSERT( box_bounds[3] == centroid_x + radius );
	TEST_ASSERT( box_bounds[4] == centroid_y + radius );
	TEST_ASSERT( box_bounds[5] == centroid_z + length/2 );

	// Check the centroid.
	Teuchos::Array<double> cylinder_centroid = cylinder.centroid();
	TEST_ASSERT( cylinder_centroid[0] == centroid_x );
	TEST_ASSERT( cylinder_centroid[1] == centroid_y );
	TEST_ASSERT( cylinder_centroid[2] == centroid_z );

	// Test some random points inside of it.
	Teuchos::Array<double> point(3);
	int num_rand = 100;
	double x_distance = 0.0;
	double y_distance = 0.0;
	double centroid_distance = 0.0;
	double tol = 1.0e-6;
	for ( int i = 0; i < num_rand; ++i )
	{
	    point[0] = 2.0 * (double) std::rand() / RAND_MAX - 1.0;
	    point[1] = 2.0 * (double) std::rand() / RAND_MAX - 1.0;
	    point[2] = 2.0 * (double) std::rand() / RAND_MAX - 1.0;

	    x_distance = centroid_x - point[0];
	    y_distance = centroid_y - point[1];
	    centroid_distance = x_distance*x_distance + y_distance*y_distance;
	    centroid_distance = std::pow( centroid_distance, 0.5 );

	    if ( centroid_distance <= radius + tol &&
		 centroid_z - length/2 <= point[2] + tol &&
		 centroid_z + length/2 >= point[2] - tol )
	    {
		TEST_ASSERT( cylinder.pointInEntity( point(), tol ) );
	    }
	    else
	    {
		TEST_ASSERT( !cylinder.pointInEntity( point(), tol ) );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Cylinder, cylinder_traits_test )
{
    using namespace DataTransferKit;

    // Make sure that PI is PI.
    double zero = 0.0;
    double pi = 2.0 * std::acos(zero);
    TEST_ASSERT( softEquivalence( pi, 3.14159, 1.0e-5 ) );

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::RCP<GeometricEntity> entity;
    for ( int i = 0; i < num_cylinders; ++i )
    {
	// Make a cylinder.
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5;
	entity = Teuchos::rcp(
	    new Cylinder( length, radius, centroid_x, centroid_y, centroid_z) );

	// Compute the measure.
	double measure = pi*radius*radius*length;
	TEST_ASSERT( softEquivalence( entity->measure(), measure, 1.0e-6 ) );

	// Compute the bounding box.
	BoundingBox box = entity->boundingBox();
	Teuchos::Tuple<double,6> box_bounds = box.getBounds();
	TEST_ASSERT( box_bounds[0] == centroid_x - radius );
	TEST_ASSERT( box_bounds[1] == centroid_y - radius );
	TEST_ASSERT( box_bounds[2] == centroid_z - length/2 );
	TEST_ASSERT( box_bounds[3] == centroid_x + radius );
	TEST_ASSERT( box_bounds[4] == centroid_y + radius );
	TEST_ASSERT( box_bounds[5] == centroid_z + length/2 );

	// Check the centroid.
	Teuchos::Array<double> cylinder_centroid = entity->centroid();
	TEST_ASSERT( cylinder_centroid[0] == centroid_x );
	TEST_ASSERT( cylinder_centroid[1] == centroid_y );
	TEST_ASSERT( cylinder_centroid[2] == centroid_z );

	// Test some random points inside of it.
	Teuchos::Array<double> point(3);
	int num_rand = 100;
	double x_distance = 0.0;
	double y_distance = 0.0;
	double centroid_distance = 0.0;
	double tol = 1.0e-6;
	for ( int i = 0; i < num_rand; ++i )
	{
	    point[0] = 2.0 * (double) std::rand() / RAND_MAX - 1.0;
	    point[1] = 2.0 * (double) std::rand() / RAND_MAX - 1.0;
	    point[2] = 2.0 * (double) std::rand() / RAND_MAX - 1.0;

	    x_distance = centroid_x - point[0];
	    y_distance = centroid_y - point[1];
	    centroid_distance = x_distance*x_distance + y_distance*y_distance;
	    centroid_distance = std::pow( centroid_distance, 0.5 );

	    if ( centroid_distance <= radius + tol &&
		 centroid_z - length/2 <= point[2] + tol &&
		 centroid_z + length/2 >= point[2] - tol )
	    {
		TEST_ASSERT( entity->pointInEntity(point(), tol) );
	    }
	    else
	    {
		TEST_ASSERT( !entity->pointInEntity(point(), tol) );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// end tstCylinder.cpp
//---------------------------------------------------------------------------//

