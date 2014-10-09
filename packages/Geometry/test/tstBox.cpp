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
#include <DTK_Entity.hpp>
#include <DTK_MappingStatus.hpp>
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

    // make a box.
    double x_min = 3.2;
    double y_min = -9.233;
    double z_min = 1.3;
    double x_max = 4.3;
    double y_max = 0.3;
    double z_max = 8.7;
    Box box(  0, 0, x_min, y_min, z_min, x_max, y_max, z_max );

    // Check Entity data.
    TEST_EQUALITY( box.name(), "DTK Box" );
    TEST_EQUALITY( box.entityType(), VOLUME );
    TEST_EQUALITY( box.id(), 0 );
    TEST_EQUALITY( box.ownerRank(), 0 );
    TEST_EQUALITY( box.physicalDimension(), 3 );
    TEST_EQUALITY( box.parametricDimension(), 3 );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds;
    box.boundingBox( box_bounds );
    TEST_ASSERT( box_bounds[0] == x_min );
    TEST_ASSERT( box_bounds[1] == y_min );
    TEST_ASSERT( box_bounds[2] == z_min );
    TEST_ASSERT( box_bounds[3] == x_max );
    TEST_ASSERT( box_bounds[4] == y_max );
    TEST_ASSERT( box_bounds[5] == z_max );

    // Compute the measure.
    TEST_FLOATING_EQUALITY( box.measure(), 77.5986, 1.0e-4 );

    // Test some random points inside of it.
    Teuchos::Array<double> point(3);
    Teuchos::Array<double> ref_point(3);
    MappingStatus status;
    Teuchos::ParameterList plist;
    plist.set<double>("Inclusion Tolerance",1.0e-12);
    bool point_inclusion = false;
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	box.safeguardMapToReferenceFrame( plist, point(), status );
	TEST_ASSERT( status.success() );

	box.mapToReferenceFrame( plist, point, ref_point(), status );
	TEST_ASSERT( status.success() );

	point_inclusion = box.checkPointInclusion(plist,ref_point());

	if ( box_bounds[0] <= ref_point[0] &&
	     box_bounds[1] <= ref_point[1] &&
	     box_bounds[2] <= ref_point[2] &&
	     box_bounds[3] >= ref_point[0] &&
	     box_bounds[4] >= ref_point[1] &&
	     box_bounds[5] >= ref_point[2] )
	{
	    TEST_ASSERT( point_inclusion );
	}
	else
	{
	    TEST_ASSERT( !point_inclusion );
	}
    }
}

//---------------------------------------------------------------------------//
// Tuple constructor test.
TEUCHOS_UNIT_TEST( Box, tuple_constructor_test )
{
    using namespace DataTransferKit;

    // make a box.
    Teuchos::Tuple<double,6> input_bounds;
    input_bounds[0] = 3.2;
    input_bounds[1] = -9.233;
    input_bounds[2] = 1.3;
    input_bounds[3] = 4.3;
    input_bounds[4] = 0.3;
    input_bounds[5] = 8.7;
    Teuchos::RCP<Entity> box = 
	Teuchos::rcp( new Box(0,0,input_bounds) );

    // Check Entity data.
    TEST_EQUALITY( box->name(), "DTK Box" );
    TEST_EQUALITY( box->entityType(), VOLUME );
    TEST_EQUALITY( box->id(), 0 );
    TEST_EQUALITY( box->ownerRank(), 0 );
    TEST_EQUALITY( box->physicalDimension(), 3 );
    TEST_EQUALITY( box->parametricDimension(), 3 );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds;
    box->boundingBox( box_bounds );
    TEST_ASSERT( box_bounds[0] == input_bounds[0] );
    TEST_ASSERT( box_bounds[1] == input_bounds[1] );
    TEST_ASSERT( box_bounds[2] == input_bounds[2] );
    TEST_ASSERT( box_bounds[3] == input_bounds[3] );
    TEST_ASSERT( box_bounds[4] == input_bounds[4] );
    TEST_ASSERT( box_bounds[5] == input_bounds[5] );

    // Compute the measure.
    TEST_FLOATING_EQUALITY( box->measure(), 77.5986, 1.0e-4 );

    // Test some random points inside of it.
    Teuchos::Array<double> point(3);
    Teuchos::Array<double> ref_point(3);
    MappingStatus status;
    Teuchos::ParameterList plist;
    plist.set<double>("Inclusion Tolerance",1.0e-12);
    bool point_inclusion = false;
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	box->safeguardMapToReferenceFrame( plist, point(), status );
	TEST_ASSERT( status.success() );

	box->mapToReferenceFrame( plist, point, ref_point(), status );
	TEST_ASSERT( status.success() );

	point_inclusion = box->checkPointInclusion(plist,ref_point());

	if ( box_bounds[0] <= ref_point[0] &&
	     box_bounds[1] <= ref_point[1] &&
	     box_bounds[2] <= ref_point[2] &&
	     box_bounds[3] >= ref_point[0] &&
	     box_bounds[4] >= ref_point[1] &&
	     box_bounds[5] >= ref_point[2] )
	{
	    TEST_ASSERT( point_inclusion );
	}
	else
	{
	    TEST_ASSERT( !point_inclusion );
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

    Box box_1( 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

    has_intersect = Box::intersectBoxes( box_1, box_2, intersection );
    TEST_ASSERT( has_intersect );
    intersection.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], 0.25 );
    TEST_EQUALITY( bounds[3], 0.75 );
    TEST_EQUALITY( bounds[4], 0.75 );
    TEST_EQUALITY( bounds[5], 0.75 );

    has_intersect = Box::intersectBoxes( box_1, box_1, intersection );
    TEST_ASSERT( has_intersect );
    intersection.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    has_intersect = Box::intersectBoxes( box_3, box_1, intersection );
    TEST_ASSERT( has_intersect );
    intersection.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 0.67 );
    TEST_EQUALITY( bounds[4], 0.67 );
    TEST_EQUALITY( bounds[5], 0.67 );

    has_intersect = Box::intersectBoxes( box_4, box_1, intersection );
    TEST_ASSERT( !has_intersect );

    has_intersect = Box::intersectBoxes( box_1, box_5, intersection );
    TEST_ASSERT( has_intersect );
    intersection.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 1.0 );
    TEST_EQUALITY( bounds[1], 1.0 );
    TEST_EQUALITY( bounds[2], 1.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    has_intersect = Box::intersectBoxes( box_5, box_1, intersection );
    TEST_ASSERT( has_intersect );
    intersection.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 1.0 );
    TEST_EQUALITY( bounds[1], 1.0 );
    TEST_EQUALITY( bounds[2], 1.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    has_intersect = Box::intersectBoxes( box_2, box_3, intersection );
    TEST_ASSERT( has_intersect );
    intersection.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], 0.25 );
    TEST_EQUALITY( bounds[3], 0.67 );
    TEST_EQUALITY( bounds[4], 0.67 );
    TEST_EQUALITY( bounds[5], 0.67 );

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

    Box box_1( 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

    Box::uniteBoxes( box_1, box_2, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    Box::uniteBoxes( box_1, box_1, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    Box::uniteBoxes( box_3, box_1, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    Box::uniteBoxes( box_4, box_1, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 1.0 );

    Box::uniteBoxes( box_1, box_5, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    Box::uniteBoxes( box_5, box_1, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    Box::uniteBoxes( box_2, box_3, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 0.75 );
    TEST_EQUALITY( bounds[4], 0.75 );
    TEST_EQUALITY( bounds[5], 0.75 );

    Box::uniteBoxes( box_2, box_4, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 0.75 );

    Box::uniteBoxes( box_2, box_5, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], 0.25 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    Box::uniteBoxes( box_3, box_5, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    Box::uniteBoxes( box_3, box_4, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 0.67 );

    Box::uniteBoxes( box_4, box_5, box_union );
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 1.0 );
    TEST_EQUALITY( bounds[1], 1.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 1.1 );
}

//---------------------------------------------------------------------------//
// Box add test.
TEUCHOS_UNIT_TEST( Box, add_test )
{
    using namespace DataTransferKit;
 
    Box box_union;
    Teuchos::Tuple<double,6> bounds;

    Box box_1( 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

    box_union = box_1 + box_2;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_1 + box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_3 + box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_4 + box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_1 + box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_5 + box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_2 + box_3;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 0.75 );
    TEST_EQUALITY( bounds[4], 0.75 );
    TEST_EQUALITY( bounds[5], 0.75 );

    box_union = box_2 + box_4;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 0.75 );

    box_union = box_2 + box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], 0.25 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_3 + box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_3 + box_4;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 0.67 );

    box_union = box_4 + box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 1.0 );
    TEST_EQUALITY( bounds[1], 1.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 1.1 );
}

//---------------------------------------------------------------------------//
// Box compound test.
TEUCHOS_UNIT_TEST( Box, compound_test )
{
    using namespace DataTransferKit;
 
    Box box_union;
    Teuchos::Tuple<double,6> bounds;

    Box box_1( 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

    box_union = box_1;
    box_union += box_2;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_1;
    box_union += box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_3;
    box_union += box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 1.0 );
    TEST_EQUALITY( bounds[4], 1.0 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_4;
    box_union += box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 1.0 );

    box_union = box_1;
    box_union += box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_5;
    box_union += box_1;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.0 );
    TEST_EQUALITY( bounds[1], 0.0 );
    TEST_EQUALITY( bounds[2], 0.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_2;
    box_union += box_3;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 0.75 );
    TEST_EQUALITY( bounds[4], 0.75 );
    TEST_EQUALITY( bounds[5], 0.75 );

    box_union = box_2;
    box_union += box_4;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 0.75 );

    box_union = box_2;
    box_union += box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 0.25 );
    TEST_EQUALITY( bounds[1], 0.25 );
    TEST_EQUALITY( bounds[2], 0.25 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_3;
    box_union += box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.0 );
    TEST_EQUALITY( bounds[3], 1.1 );
    TEST_EQUALITY( bounds[4], 1.1 );
    TEST_EQUALITY( bounds[5], 1.1 );

    box_union = box_3;
    box_union += box_4;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], -1.0 );
    TEST_EQUALITY( bounds[1], -1.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 0.67 );

    box_union = box_4;
    box_union += box_5;
    box_union.boundingBox( bounds );
    TEST_EQUALITY( bounds[0], 1.0 );
    TEST_EQUALITY( bounds[1], 1.0 );
    TEST_EQUALITY( bounds[2], -1.2 );
    TEST_EQUALITY( bounds[3], 5.6 );
    TEST_EQUALITY( bounds[4], 7.8 );
    TEST_EQUALITY( bounds[5], 1.1 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Box, inclusion_test )
{
    using namespace DataTransferKit;

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::RCP<Entity> box;
    for ( int i = 0; i < num_boxes; ++i )
    {
	// Make a box.
	double x_min = -(double) std::rand() / RAND_MAX;
	double y_min = -(double) std::rand() / RAND_MAX;
	double z_min = -(double) std::rand() / RAND_MAX;
	double x_max =  (double) std::rand() / RAND_MAX;
	double y_max =  (double) std::rand() / RAND_MAX;
	double z_max =  (double) std::rand() / RAND_MAX;
	box = Teuchos::rcp(
	    new Box( 1, 1, x_min, y_min, z_min, x_max, y_max, z_max ) );

	// Compute the measure.
	double measure = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	TEST_FLOATING_EQUALITY( box->measure(), measure, 1.0e-6 );

	// Check the centroid.
	Teuchos::ArrayView<const double> box_centroid; 
	box->centroid( box_centroid );
	TEST_ASSERT( box_centroid[0] == (x_max+x_min)/2.0 );
	TEST_ASSERT( box_centroid[1] == (y_max+y_min)/2.0 );
	TEST_ASSERT( box_centroid[2] == (z_max+z_min)/2.0 );

	// Test some random points inside of it.
	int num_rand = 100;
	double tol = 1.0e-6;
	Teuchos::Array<double> point(3);
	Teuchos::Array<double> ref_point(3);
	MappingStatus status;
	Teuchos::ParameterList plist;
	plist.set<double>("Inclusion Tolerance",tol);
	bool point_inclusion = false;
	for ( int i = 0; i < num_rand; ++i )
	{
	    point[0] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[1] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;
	    point[2] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;

	    box->safeguardMapToReferenceFrame( plist, point(), status );
	    TEST_ASSERT( status.success() );

	    box->mapToReferenceFrame( plist, point, ref_point(), status );
	    TEST_ASSERT( status.success() );

	    point_inclusion = box->checkPointInclusion(plist,ref_point());

	    if ( x_min - tol <= point[0] && point[0] <= x_max + tol &&
		 y_min - tol <= point[1] && point[1] <= y_max + tol &&
		 z_min - tol <= point[2] && point[2] <= z_max + tol )
	    {
		TEST_ASSERT( point_inclusion );
	    }
	    else
	    {
		TEST_ASSERT( !point_inclusion );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Box, communication_test )
{
    using namespace DataTransferKit;

    // Register the box class to use the abstract compile-time interfaces.
    AbstractObjectRegistry<Entity,Box>::registerDerivedClasses();

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_default = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm_default->getRank();

    // make a box.
    double x_min = 3.2;
    double y_min = -9.233;
    double z_min = 1.3;
    double x_max = 4.3;
    double y_max = 0.3;
    double z_max = 8.7;

    Box box;
    Entity entity;
    if ( 0 == comm_rank )
    {
	box = Box(  0, 0, x_min, y_min, z_min, x_max, y_max, z_max );
	entity = Box(0, 0, x_min, y_min, z_min, x_max, y_max, z_max);
    }

    // Directly serialize the subclass.
    Teuchos::broadcast( *comm_default, 0, Teuchos::Ptr<Box>(&box) );

    // Broadcast the box with indirect serialization through the geometric
    // entity api.
    Teuchos::broadcast( *comm_default, 0, Teuchos::Ptr<Entity>(&entity) );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds;
    box.boundingBox( box_bounds );
    TEST_EQUALITY( box_bounds[0], x_min );
    TEST_EQUALITY( box_bounds[1], y_min );
    TEST_EQUALITY( box_bounds[2], z_min );
    TEST_EQUALITY( box_bounds[3], x_max );
    TEST_EQUALITY( box_bounds[4], y_max );
    TEST_EQUALITY( box_bounds[5], z_max );

    Teuchos::Tuple<double,6> entity_bounds;
    entity.boundingBox( entity_bounds );
    TEST_EQUALITY( entity_bounds[0], x_min );
    TEST_EQUALITY( entity_bounds[1], y_min );
    TEST_EQUALITY( entity_bounds[2], z_min );
    TEST_EQUALITY( entity_bounds[3], x_max );
    TEST_EQUALITY( entity_bounds[4], y_max );
    TEST_EQUALITY( entity_bounds[5], z_max );

    // Compute the measure.
    TEST_FLOATING_EQUALITY( box.measure(), 77.5986, 1.0e-4 );
    TEST_FLOATING_EQUALITY( entity.measure(), 77.5986, 1.0e-4 );

    // Test some random points inside of the box.
    Teuchos::Array<double> point(3);
    Teuchos::Array<double> ref_point(3);
    MappingStatus status;
    Teuchos::ParameterList plist;
    plist.set<double>("Inclusion Tolerance",1.0e-12);
    bool point_inclusion = false;
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	box.safeguardMapToReferenceFrame( plist, point(), status );
	TEST_ASSERT( status.success() );

	box.mapToReferenceFrame( plist, point, ref_point(), status );
	TEST_ASSERT( status.success() );

	point_inclusion = box.checkPointInclusion(plist,ref_point());

	if ( box_bounds[0] <= ref_point[0] &&
	     box_bounds[1] <= ref_point[1] &&
	     box_bounds[2] <= ref_point[2] &&
	     box_bounds[3] >= ref_point[0] &&
	     box_bounds[4] >= ref_point[1] &&
	     box_bounds[5] >= ref_point[2] )
	{
	    TEST_ASSERT( point_inclusion );
	}
	else
	{
	    TEST_ASSERT( !point_inclusion );
	}
    }

    // Test some random points inside of the entity.
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	entity.safeguardMapToReferenceFrame( plist, point(), status );
	TEST_ASSERT( status.success() );

	entity.mapToReferenceFrame( plist, point, ref_point(), status );
	TEST_ASSERT( status.success() );

	point_inclusion = entity.checkPointInclusion(plist,ref_point());

	if ( entity_bounds[0] <= ref_point[0] &&
	     entity_bounds[1] <= ref_point[1] &&
	     entity_bounds[2] <= ref_point[2] &&
	     entity_bounds[3] >= ref_point[0] &&
	     entity_bounds[4] >= ref_point[1] &&
	     entity_bounds[5] >= ref_point[2] )
	{
	    TEST_ASSERT( point_inclusion );
	}
	else
	{
	    TEST_ASSERT( !point_inclusion );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstBox.cpp
//---------------------------------------------------------------------------//

