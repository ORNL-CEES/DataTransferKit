//---------------------------------------------------------------------------//
/*!
 * \file tstBoundingBox.cpp
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
// Helper Function
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
// Global test variables.
//---------------------------------------------------------------------------//
int num_rand = 1000;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Default constructor test.
TEUCHOS_UNIT_TEST( BoundingBox, default_constructor_test )
{
    using namespace DataTransferKit;

    // Make a bounding box.
    BoundingBox box( 3.2, -9.233, 1.3, 4.3, 0.3, 8.7 );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds = box.getBounds();
    TEST_ASSERT( box_bounds[0] == 3.2 );
    TEST_ASSERT( box_bounds[1] == -9.233 );
    TEST_ASSERT( box_bounds[2] == 1.3 );
    TEST_ASSERT( box_bounds[3] == 4.3 );
    TEST_ASSERT( box_bounds[4] == 0.3 );
    TEST_ASSERT( box_bounds[5] == 8.7 );

    // Compute the volume.
    TEST_ASSERT( softEquivalence( box.volume( 3 ), 77.5986, 1.0e-4 ) );

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
TEUCHOS_UNIT_TEST( BoundingBox, tuple_constructor_test )
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
    BoundingBox box( input_bounds );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds = box.getBounds();
    TEST_ASSERT( box_bounds[0] == 3.2 );
    TEST_ASSERT( box_bounds[1] == -9.233 );
    TEST_ASSERT( box_bounds[2] == 1.3 );
    TEST_ASSERT( box_bounds[3] == 4.3 );
    TEST_ASSERT( box_bounds[4] == 0.3 );
    TEST_ASSERT( box_bounds[5] == 8.7 );

    // Compute the volume.
    TEST_ASSERT( softEquivalence( box.volume( 3 ), 77.5986, 1.0e-4 ) );

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
TEUCHOS_UNIT_TEST( BoundingBox, intersection_test )
{
    using namespace DataTransferKit;
 
    bool has_intersect;
    BoundingBox intersection;
    Teuchos::Tuple<double,6> bounds;

    BoundingBox box_1( 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    BoundingBox box_2( 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    BoundingBox box_3( -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    BoundingBox box_4( 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    BoundingBox box_5( 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

    has_intersect = BoundingBox::intersectBoxes( box_1, box_2, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.25 );
    TEST_ASSERT( bounds[1] == 0.25 );
    TEST_ASSERT( bounds[2] == 0.25 );
    TEST_ASSERT( bounds[3] == 0.75 );
    TEST_ASSERT( bounds[4] == 0.75 );
    TEST_ASSERT( bounds[5] == 0.75 );

    has_intersect = BoundingBox::intersectBoxes( box_1, box_1, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    has_intersect = BoundingBox::intersectBoxes( box_3, box_1, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.0 );
    TEST_ASSERT( bounds[1] == 0.0 );
    TEST_ASSERT( bounds[2] == 0.0 );
    TEST_ASSERT( bounds[3] == 0.67 );
    TEST_ASSERT( bounds[4] == 0.67 );
    TEST_ASSERT( bounds[5] == 0.67 );

    has_intersect = BoundingBox::intersectBoxes( box_4, box_1, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = BoundingBox::intersectBoxes( box_1, box_5, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 1.0 );
    TEST_ASSERT( bounds[1] == 1.0 );
    TEST_ASSERT( bounds[2] == 1.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    has_intersect = BoundingBox::intersectBoxes( box_5, box_1, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 1.0 );
    TEST_ASSERT( bounds[1] == 1.0 );
    TEST_ASSERT( bounds[2] == 1.0 );
    TEST_ASSERT( bounds[3] == 1.0 );
    TEST_ASSERT( bounds[4] == 1.0 );
    TEST_ASSERT( bounds[5] == 1.0 );

    has_intersect = BoundingBox::intersectBoxes( box_2, box_3, intersection );
    TEST_ASSERT( has_intersect );
    bounds = intersection.getBounds();
    TEST_ASSERT( bounds[0] == 0.25 );
    TEST_ASSERT( bounds[1] == 0.25 );
    TEST_ASSERT( bounds[2] == 0.25 );
    TEST_ASSERT( bounds[3] == 0.67 );
    TEST_ASSERT( bounds[4] == 0.67 );
    TEST_ASSERT( bounds[5] == 0.67 );

    has_intersect = BoundingBox::intersectBoxes( box_2, box_4, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = BoundingBox::intersectBoxes( box_2, box_5, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = BoundingBox::intersectBoxes( box_3, box_5, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = BoundingBox::intersectBoxes( box_3, box_4, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = BoundingBox::intersectBoxes( box_4, box_5, intersection );
    TEST_ASSERT( !has_intersect );
}

//---------------------------------------------------------------------------//
// end tstBoundingBox.cpp
//---------------------------------------------------------------------------//

