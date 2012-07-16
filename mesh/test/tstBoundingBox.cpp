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

    // Test some points inside of it.
    double point_0[3] = { 3.7, -4, 5.4 };
    double point_1[3] = { 4.25, -7.99, 8.3 };
    double point_2[3] = { 5.4, -3, 9.4 };
    double point_3[3] = { 2.7, 0.4, 8.3 };
    double point_4[3] = { 3.2, -9.233, 1.3 };
    double point_5[3] = { 3.2, 0.3, 1.3 };
    TEST_ASSERT( box.pointInBox( point_0 ) );
    TEST_ASSERT( box.pointInBox( point_1 ) );
    TEST_ASSERT( !box.pointInBox( point_2 ) );
    TEST_ASSERT( !box.pointInBox( point_3 ) );
    TEST_ASSERT( box.pointInBox( point_4 ) );
    TEST_ASSERT( box.pointInBox( point_5 ) );
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

    // Test some points inside of it.
    double point_0[3] = { 3.7, -4, 5.4 };
    double point_1[3] = { 4.25, -7.99, 8.3 };
    double point_2[3] = { 5.4, -3, 9.4 };
    double point_3[3] = { 2.7, 0.4, 8.3 };
    double point_4[3] = { 3.2, -9.233, 1.3 };
    double point_5[3] = { 3.2, 0.3, 1.3 };
    TEST_ASSERT( box.pointInBox( point_0 ) );
    TEST_ASSERT( box.pointInBox( point_1 ) );
    TEST_ASSERT( !box.pointInBox( point_2 ) );
    TEST_ASSERT( !box.pointInBox( point_3 ) );
    TEST_ASSERT( box.pointInBox( point_4 ) );
    TEST_ASSERT( box.pointInBox( point_5 ) );
}

//---------------------------------------------------------------------------//
// Serialization test.
TEUCHOS_UNIT_TEST( BoundingBox, serialization_test )
{
    using namespace DataTransferKit;

    // Comm setup.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Make a bounding box on each process.
    Teuchos::Array<BoundingBox> boxes( my_size );
    boxes[my_rank] = BoundingBox( my_rank, my_rank, my_rank,
				  my_rank+1.0, my_rank+1.0, my_rank+1.0 );

    // Reduce the boxes to an array on all processes.
    Teuchos::reduceAll<int,BoundingBox>( *comm,
					 Teuchos::REDUCE_SUM,
					 boxes.size(),
					 &boxes[0],
					 &boxes[0] );

    // Check the reduction with rank dependent coordinates.
    double rank = 0.0;
    Teuchos::Array<BoundingBox>::const_iterator box_iterator;
    for ( box_iterator = boxes.begin(); 
	  box_iterator != boxes.end(); 
	  ++box_iterator )
    {
	double point[3] = { rank+0.5, rank+0.5, rank+0.5 };
	TEST_ASSERT( box_iterator->pointInBox( point ) );

	rank += 1.0;
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

