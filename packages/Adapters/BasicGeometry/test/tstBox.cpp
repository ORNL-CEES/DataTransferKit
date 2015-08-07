//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
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
    Box box(  0, 0, 0, x_min, y_min, z_min, x_max, y_max, z_max );
    BasicGeometryEntity box_entity = box;
    Entity entity = box;

    // Check Entity data.
    TEST_EQUALITY( box.id(), 0 );
    TEST_EQUALITY( box.ownerRank(), 0 );
    TEST_ASSERT( box.inBlock(0) );
    TEST_ASSERT( !box.onBoundary(0) );
    TEST_EQUALITY( box.topologicalDimension(), 3 );
    TEST_EQUALITY( box.physicalDimension(), 3 );

    TEST_EQUALITY( box_entity.id(), 0 );
    TEST_EQUALITY( box_entity.ownerRank(), 0 );
    TEST_ASSERT( box_entity.inBlock(0) );
    TEST_ASSERT( !box_entity.onBoundary(0) );
    TEST_EQUALITY( box_entity.topologicalDimension(), 3 );
    TEST_EQUALITY( box_entity.physicalDimension(), 3 );

    TEST_EQUALITY( entity.id(), 0 );
    TEST_EQUALITY( entity.ownerRank(), 0 );
    TEST_ASSERT( entity.inBlock(0) );
    TEST_ASSERT( !entity.onBoundary(0) );
    TEST_EQUALITY( entity.topologicalDimension(), 3 );
    TEST_EQUALITY( entity.physicalDimension(), 3 );

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
    Teuchos::ParameterList plist;
    double tol = 1.0e-12;
    bool point_inclusion = false;
    bool map_ok = false;
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	// Box API
	map_ok = box.mapToReferenceFrame( point, ref_point() );
	TEST_ASSERT( map_ok );
	point_inclusion = box.checkPointInclusion(tol,ref_point());
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

	// BasicGeometryEntity API
	map_ok = box_entity.mapToReferenceFrame( point, ref_point() );
	TEST_ASSERT( map_ok );
	point_inclusion = box_entity.checkPointInclusion(tol,ref_point());
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
    Teuchos::RCP<Box> box = 
	Teuchos::rcp( new Box(0,0,0,input_bounds) );
    Teuchos::RCP<Entity> entity = box;

    // Check Entity data.
    TEST_EQUALITY( entity->id(), 0 );
    TEST_EQUALITY( entity->ownerRank(), 0 );
    TEST_ASSERT( entity->inBlock(0) );
    TEST_ASSERT( !entity->onBoundary(0) );
    TEST_EQUALITY( entity->topologicalDimension(), 3 );
    TEST_EQUALITY( entity->physicalDimension(), 3 );

    // Check the bounds.
    Teuchos::Tuple<double,6> box_bounds;
    entity->boundingBox( box_bounds );
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
    Teuchos::ParameterList plist;
    double tol = 1.0e-12;
    bool point_inclusion = false;
    bool map_ok = false;
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	map_ok = box->mapToReferenceFrame( point, ref_point() );
	TEST_ASSERT( map_ok );
	point_inclusion = box->checkPointInclusion(tol,ref_point());

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

    Box box_1( 0, 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

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

    Box box_1( 0, 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

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

    Box box_1( 0, 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

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

    Box box_1( 0, 0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    Box box_2( 0, 0, 0, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
    Box box_3( 0, 0, 0, -1.0, -1.0, -1.0, 0.67, 0.67, 0.67);
    Box box_4( 0, 0, 0, 4.3, 6.2, -1.2, 5.6, 7.8, -0.8 );
    Box box_5( 0, 0, 0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1 );

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
// end tstBox.cpp
//---------------------------------------------------------------------------//

