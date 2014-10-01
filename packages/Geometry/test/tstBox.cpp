//---------------------------------------------------------------------------//
/*!
 * \file tstBox.cpp
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

#include <DTK_Box.hpp>
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
// Global test variables.
//---------------------------------------------------------------------------//
int num_rand = 1000;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Default constructor test.
TEUCHOS_UNIT_TEST( Box, default_constructor_test )
{
    using namespace DataTransferKit;

    // Make a bounding box.
    Box box( 3.2, -9.233, 1.3, 4.3, 0.3, 8.7 );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds = box.getBounds();
    TEST_ASSERT( box_bounds[0] == 3.2 );
    TEST_ASSERT( box_bounds[1] == -9.233 );
    TEST_ASSERT( box_bounds[2] == 1.3 );
    TEST_ASSERT( box_bounds[3] == 4.3 );
    TEST_ASSERT( box_bounds[4] == 0.3 );
    TEST_ASSERT( box_bounds[5] == 8.7 );

    // Compute the volume.
    TEST_FLOATING_EQUALITY( box.volume( 3 ), 77.5986, 1.0e-4 );

    // Test some random points inside of it.
    Teuchos::Array<double> point(3);
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	if ( box_bounds[0] <= point[0] &&
	     box_bounds[1] <= point[1] &&
	     box_bounds[2] <= point[2] &&
	     box_bounds[3] >= point[0] &&
	     box_bounds[4] >= point[1] &&
	     box_bounds[5] >= point[2] )
	{
	    TEST_ASSERT( box.pointInBox( point ) );
	}
	else
	{
	    TEST_ASSERT( !box.pointInBox( point ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Tuple constructor test.
TEUCHOS_UNIT_TEST( Box, tuple_constructor_test )
{
    using namespace DataTransferKit;

    // Make a bounding box.
    Teuchos::Tuple<double,6> input_bounds;
    input_bounds[0] = 3.2;
    input_bounds[1] = -9.233;
    input_bounds[2] = 1.3;
    input_bounds[3] = 4.3;
    input_bounds[4] = 0.3;
    input_bounds[5] = 8.7;
    Box box( input_bounds );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds = box.getBounds();
    TEST_ASSERT( box_bounds[0] == 3.2 );
    TEST_ASSERT( box_bounds[1] == -9.233 );
    TEST_ASSERT( box_bounds[2] == 1.3 );
    TEST_ASSERT( box_bounds[3] == 4.3 );
    TEST_ASSERT( box_bounds[4] == 0.3 );
    TEST_ASSERT( box_bounds[5] == 8.7 );

    // Compute the volume.
    TEST_FLOATING_EQUALITY( box.volume( 3 ), 77.5986, 1.0e-4 );

    // Test some random points inside of it.
    Teuchos::Array<double> point(3);
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	if ( box_bounds[0] <= point[0] &&
	     box_bounds[1] <= point[1] &&
	     box_bounds[2] <= point[2] &&
	     box_bounds[3] >= point[0] &&
	     box_bounds[4] >= point[1] &&
	     box_bounds[5] >= point[2] )
	{
	    TEST_ASSERT( box.pointInBox( point ) );
	}
	else
	{
	    TEST_ASSERT( !box.pointInBox( point ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Box intersection test.
TEUCHOS_UNIT_TEST( Box, intersection_test )
{
    using namespace DataTransferKit;
 
    bool has_intersect;
    Box intersection;
    Teuchos::Tuple<double,6> bounds;

    Box box_1( 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

    has_intersect = Box::intersectBoxes( box_1, box_2, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.25 );
    TEST_ASSERT( bounds[1] == 0.25 );
    TEST_ASSERT( bounds[2] == 0.25 );
    TEST_ASSERT( bounds[3] == 0.75 );
    TEST_ASSERT( bounds[4] == 0.75 );
    TEST_ASSERT( bounds[5] == 0.75 );

    has_intersect = Box::intersectBoxes( box_1, box_1, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    has_intersect = Box::intersectBoxes( box_3, box_1, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 0.67 );
    TEST_ASSERT( bounds[4] == 0.67 );
    TEST_ASSERT( bounds[5] == 0.67 );

    has_intersect = Box::intersectBoxes( box_4, box_1, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = Box::intersectBoxes( box_1, box_5, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 1.0 );
    TEST_ASSERT( bounds[1] == 1.0 );
    TEST_ASSERT( bounds[2] == 1.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    has_intersect = Box::intersectBoxes( box_5, box_1, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 1.0 );
    TEST_ASSERT( bounds[1] == 1.0 );
    TEST_ASSERT( bounds[2] == 1.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    has_intersect = Box::intersectBoxes( box_2, box_3, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.25 );
    TEST_ASSERT( bounds[1] == 0.25 );
    TEST_ASSERT( bounds[2] == 0.25 );
    TEST_ASSERT( bounds[3] == 0.67 );
    TEST_ASSERT( bounds[4] == 0.67 );
    TEST_ASSERT( bounds[5] == 0.67 );

    has_intersect = Box::intersectBoxes( box_2, box_4, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = Box::intersectBoxes( box_2, box_5, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = Box::intersectBoxes( box_3, box_5, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = Box::intersectBoxes( box_3, box_4, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = Box::intersectBoxes( box_4, box_5, intersection );
    TEST_ASSERT( !has_intersect );
}

//---------------------------------------------------------------------------//
// Box union test.
TEUCHOS_UNIT_TEST( Box, union_test )
{
    using namespace DataTransferKit;
 
    Box box_union;
    Teuchos::Tuple<double,6> bounds;

    Box box_1( 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

    Box::uniteBoxes( box_1, box_2, box_union );
    bounds = box_union.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    Box::uniteBoxes( box_1, box_1, box_union );
    bounds = box_union.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    Box::uniteBoxes( box_3, box_1, box_union );
    bounds = box_union.getBounds();
    TEST_ASSERT( bounds[0] == -1.0 );
    TEST_ASSERT( bounds[1] == -1.0 );
    TEST_ASSERT( bounds[2] == -1.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    Box::uniteBoxes( box_4, box_1, box_union );
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == -1.2 );
    TEST_ASSERT( bounds[3] == 5.6 );
    TEST_ASSERT( bounds[4] == 7.8 );
    TEST_ASSERT( bounds[5] == 1.0 );

    Box::uniteBoxes( box_1, box_5, box_union );
    bounds = box_union.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 1.1 );
    TEST_ASSERT( bounds[4] == 1.1 );
    TEST_ASSERT( bounds[5] == 1.1 );

    Box::uniteBoxes( box_5, box_1, box_union );
    bounds = box_union.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 1.1 );
    TEST_ASSERT( bounds[4] == 1.1 );
    TEST_ASSERT( bounds[5] == 1.1 );

    Box::uniteBoxes( box_2, box_3, box_union );
    bounds = box_union.getBounds();
    TEST_ASSERT( bounds[0] == -1.0 );
    TEST_ASSERT( bounds[1] == -1.0 );
    TEST_ASSERT( bounds[2] == -1.0 );
    TEST_ASSERT( bounds[3] == 0.75 );
    TEST_ASSERT( bounds[4] == 0.75 );
    TEST_ASSERT( bounds[5] == 0.75 );

    Box::uniteBoxes( box_2, box_4, box_union );
    TEST_ASSERT( bounds[0] == 0.25 );
    TEST_ASSERT( bounds[1] == 0.25 );
    TEST_ASSERT( bounds[2] == -1.2 );
    TEST_ASSERT( bounds[3] == 5.6 );
    TEST_ASSERT( bounds[4] == 7.8 );
    TEST_ASSERT( bounds[5] == 0.75 );

    Box::uniteBoxes( box_2, box_5, box_union );
    TEST_ASSERT( bounds[0] == 0.25 );
    TEST_ASSERT( bounds[1] == 0.25 );
    TEST_ASSERT( bounds[2] == 0.25 );
    TEST_ASSERT( bounds[3] == 1.1 );
    TEST_ASSERT( bounds[4] == 1.1 );
    TEST_ASSERT( bounds[5] == 1.1 );

    Box::uniteBoxes( box_3, box_5, box_union );
    TEST_ASSERT( bounds[0] == -1.0 );
    TEST_ASSERT( bounds[1] == -1.0 );
    TEST_ASSERT( bounds[2] == -1.2 );
    TEST_ASSERT( bounds[3] == 1.1 );
    TEST_ASSERT( bounds[4] == 1.1 );
    TEST_ASSERT( bounds[5] == 1.1 );

    Box::uniteBoxes( box_3, box_4, box_union );
    TEST_ASSERT( bounds[0] == -1.0 );
    TEST_ASSERT( bounds[1] == -1.0 );
    TEST_ASSERT( bounds[2] == -1.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    Box::uniteBoxes( box_4, box_5, box_union );
    TEST_ASSERT( bounds[0] == 1.0 );
    TEST_ASSERT( bounds[1] == 1.0 );
    TEST_ASSERT( bounds[2] == -1.2 );
    TEST_ASSERT( bounds[3] == 5.6 );
    TEST_ASSERT( bounds[4] == 7.8 );
    TEST_ASSERT( bounds[5] == 1.1 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Box, box_test )
{
    using namespace DataTransferKit;

    // Build a series of random boxes.
    int num_boxes = 100;
    for ( int i = 0; i < num_boxes; ++i )
    {
	// Make a box.
	double x_min = -(double) std::rand() / RAND_MAX;
	double y_min = -(double) std::rand() / RAND_MAX;
	double z_min = -(double) std::rand() / RAND_MAX;
	double x_max =  (double) std::rand() / RAND_MAX;
	double y_max =  (double) std::rand() / RAND_MAX;
	double z_max =  (double) std::rand() / RAND_MAX;
	Box box( 1, 1, x_min, y_min, z_min, x_max, y_max, z_max );

	// Check the geometric bounds.
	Teuchos::Tuple<double,6> box_bounds = box.getBounds();
	TEST_ASSERT( box_bounds[0] == x_min );
	TEST_ASSERT( box_bounds[1] == y_min );
	TEST_ASSERT( box_bounds[2] == z_min );
	TEST_ASSERT( box_bounds[3] == x_max );
	TEST_ASSERT( box_bounds[4] == y_max );
	TEST_ASSERT( box_bounds[5] == z_max );

	// Compute the measure.
	double measure = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	TEST_FLOATING_EQUALITY( box.measure(), measure, 1.0e-6 );

	// Compute the bounding box.
	Box bounding_box = box.boundingBox();
	Teuchos::Tuple<double,6> bounding_box_bounds = bounding_box.getBounds();
	TEST_ASSERT( bounding_box_bounds[0] == x_min );
	TEST_ASSERT( bounding_box_bounds[1] == y_min );
	TEST_ASSERT( bounding_box_bounds[2] == z_min );
	TEST_ASSERT( bounding_box_bounds[3] == x_max );
	TEST_ASSERT( bounding_box_bounds[4] == y_max );
	TEST_ASSERT( bounding_box_bounds[5] == z_max );

	// Check the centroid.
	Teuchos::Array<double> box_centroid = box.centroid();
	TEST_ASSERT( box_centroid[0] == (x_max+x_min)/2.0 );
	TEST_ASSERT( box_centroid[1] == (y_max+y_min)/2.0 );
	TEST_ASSERT( box_centroid[2] == (z_max+z_min)/2.0 );

	// Test some random points inside of it.
	Teuchos::Array<double> point(3);
	int num_rand = 100;
	double tol = 1.0e-6;
	for ( int i = 0; i < num_rand; ++i )
	{
	    point[0] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[1] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[2] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;

	    if ( x_min - tol <= point[0] && point[0] <= x_max + tol &&
		 y_min - tol <= point[1] && point[1] <= y_max + tol &&
		 z_min - tol <= point[2] && point[2] <= z_max + tol )
	    {
		TEST_ASSERT( box.pointInEntity( point(), tol ) );
	    }
	    else
	    {
		TEST_ASSERT( !box.pointInEntity( point(), tol ) );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Box, box_traits_test )
{
    using namespace DataTransferKit;

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::RCP<GeometricEntity> entity;
    for ( int i = 0; i < num_boxes; ++i )
    {
	// Make a box.
	double x_min = -(double) std::rand() / RAND_MAX;
	double y_min = -(double) std::rand() / RAND_MAX;
	double z_min = -(double) std::rand() / RAND_MAX;
	double x_max =  (double) std::rand() / RAND_MAX;
	double y_max =  (double) std::rand() / RAND_MAX;
	double z_max =  (double) std::rand() / RAND_MAX;
	entity = Teuchos::rcp(
	    new Box( x_min, y_min, z_min, x_max, y_max, z_max ) );

	// Compute the measure.
	double measure = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	TEST_FLOATING_EQUALITY( entity->measure(), measure, 1.0e-6 );

	// Compute the bounding box.
	Box bounding_box = entity->boundingBox();
	Teuchos::Tuple<double,6> bounding_box_bounds = bounding_box.getBounds();
	TEST_ASSERT( bounding_box_bounds[0] == x_min );
	TEST_ASSERT( bounding_box_bounds[1] == y_min );
	TEST_ASSERT( bounding_box_bounds[2] == z_min );
	TEST_ASSERT( bounding_box_bounds[3] == x_max );
	TEST_ASSERT( bounding_box_bounds[4] == y_max );
	TEST_ASSERT( bounding_box_bounds[5] == z_max );

	// Check the centroid.
	Teuchos::Array<double> box_centroid = entity->centroid();
	TEST_ASSERT( box_centroid[0] == (x_max+x_min)/2.0 );
	TEST_ASSERT( box_centroid[1] == (y_max+y_min)/2.0 );
	TEST_ASSERT( box_centroid[2] == (z_max+z_min)/2.0 );

	// Test some random points inside of it.
	Teuchos::Array<double> point(3);
	int num_rand = 100;
	double tol = 1.0e-6;
	for ( int i = 0; i < num_rand; ++i )
	{
	    point[0] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[1] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[2] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;

	    if ( x_min - tol <= point[0] && point[0] <= x_max + tol &&
		 y_min - tol <= point[1] && point[1] <= y_max + tol &&
		 z_min - tol <= point[2] && point[2] <= z_max + tol )
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
// end tstBox.cpp
//---------------------------------------------------------------------------//

