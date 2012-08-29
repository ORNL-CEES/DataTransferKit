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
#include <DTK_Box.hpp>
#include <DTK_GeometryTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
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
	Box box( x_min, y_min, z_min, x_max, y_max, z_max );

	// Check the geometric bounds.
	Teuchos::Tuple<double,6> box_bounds = box.getBounds();
	TEST_ASSERT( box_bounds[0] == x_min );
	TEST_ASSERT( box_bounds[1] == y_min );
	TEST_ASSERT( box_bounds[2] == z_min );
	TEST_ASSERT( box_bounds[3] == x_max );
	TEST_ASSERT( box_bounds[4] == y_max );
	TEST_ASSERT( box_bounds[5] == z_max );

	// Compute the volume.
	double volume = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	TEST_ASSERT( softEquivalence( box.volume(), volume, 1.0e-6 ) );

	// Compute the bounding box.
	BoundingBox bounding_box = box.boundingBox();
	Teuchos::Tuple<double,6> bounding_box_bounds = bounding_box.getBounds();
	TEST_ASSERT( bounding_box_bounds[0] == x_min );
	TEST_ASSERT( bounding_box_bounds[1] == y_min );
	TEST_ASSERT( bounding_box_bounds[2] == z_min );
	TEST_ASSERT( bounding_box_bounds[3] == x_max );
	TEST_ASSERT( bounding_box_bounds[4] == y_max );
	TEST_ASSERT( bounding_box_bounds[5] == z_max );

	// Test some random points inside of it.
	Teuchos::Array<double> point(3);
	int num_rand = 100;
	for ( int i = 0; i < num_rand; ++i )
	{
	    point[0] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[1] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[2] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;

	    if ( x_min <= point[0] && point[0] <= x_max &&
		 y_min <= point[1] && point[1] <= y_max &&
		 z_min <= point[2] && point[2] <= z_max )
	    {
		TEST_ASSERT( box.pointInBox( point ) );
	    }
	    else
	    {
		TEST_ASSERT( !box.pointInBox( point ) );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Box, box_traits_test )
{
    using namespace DataTransferKit;
    typedef GeometryTraits<Box> GT;

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
	Box box( x_min, y_min, z_min, x_max, y_max, z_max );

	// Compute the volume.
	double volume = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	TEST_ASSERT( softEquivalence( GT::measure(box), volume, 1.0e-6 ) );

	// Compute the bounding box.
	BoundingBox bounding_box = GT::boundingBox(box);
	Teuchos::Tuple<double,6> bounding_box_bounds = bounding_box.getBounds();
	TEST_ASSERT( bounding_box_bounds[0] == x_min );
	TEST_ASSERT( bounding_box_bounds[1] == y_min );
	TEST_ASSERT( bounding_box_bounds[2] == z_min );
	TEST_ASSERT( bounding_box_bounds[3] == x_max );
	TEST_ASSERT( bounding_box_bounds[4] == y_max );
	TEST_ASSERT( bounding_box_bounds[5] == z_max );

	// Test some random points inside of it.
	Teuchos::Array<double> point(3);
	int num_rand = 100;
	for ( int i = 0; i < num_rand; ++i )
	{
	    point[0] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[1] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[2] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;

	    if ( x_min <= point[0] && point[0] <= x_max &&
		 y_min <= point[1] && point[1] <= y_max &&
		 z_min <= point[2] && point[2] <= z_max )
	    {
		TEST_ASSERT( GT::pointInGeometry( box, point ) );
	    }
	    else
	    {
		TEST_ASSERT( !GT::pointInGeometry( box, point ) );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// end tstBox.cpp
//---------------------------------------------------------------------------//

