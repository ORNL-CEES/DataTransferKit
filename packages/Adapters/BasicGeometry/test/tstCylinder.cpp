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

#include <DTK_Cylinder.hpp>
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

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Cylinder, cylinder_test )
{
    using namespace DataTransferKit;

    // Make sure that PI is PI.
    double zero = 0.0;
    double pi = 2.0 * std::acos(zero);
    TEST_FLOATING_EQUALITY( pi, 3.14159, 1.0e-5 );

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
	Cylinder cylinder( 
	    i, i, i, length, radius, centroid_x, centroid_y, centroid_z );

	// Check the cylinder.
	TEST_EQUALITY( cylinder.topologicalDimension(), 3 );
	TEST_EQUALITY( Teuchos::as<int>(cylinder.id()), i );
	TEST_EQUALITY( cylinder.ownerRank(), i );
	TEST_ASSERT( cylinder.inBlock(i) );
	TEST_ASSERT( !cylinder.inBlock(i+1) );
	TEST_ASSERT( !cylinder.onBoundary(i) );
	TEST_EQUALITY( cylinder.length(), length );
	TEST_EQUALITY( cylinder.radius(), radius );

	// Compute the measure.
	double measure = pi*radius*radius*length;
	TEST_FLOATING_EQUALITY( cylinder.measure(), measure, 1.0e-6 );

	// Compute the bounding box.
	Teuchos::Tuple<double,6> box_bounds;
	cylinder.boundingBox( box_bounds );
	TEST_EQUALITY( box_bounds[0], centroid_x - radius );
	TEST_EQUALITY( box_bounds[1], centroid_y - radius );
	TEST_EQUALITY( box_bounds[2], centroid_z - length/2 );
	TEST_EQUALITY( box_bounds[3], centroid_x + radius );
	TEST_EQUALITY( box_bounds[4], centroid_y + radius );
	TEST_EQUALITY( box_bounds[5], centroid_z + length/2 );

	// Check the centroid.
	Teuchos::Array<double> cylinder_centroid(3);
	cylinder.centroid( cylinder_centroid() );
	TEST_EQUALITY( cylinder_centroid[0], centroid_x );
	TEST_EQUALITY( cylinder_centroid[1], centroid_y );
	TEST_EQUALITY( cylinder_centroid[2], centroid_z );

	// Test some random points inside of it.
	Teuchos::Array<double> point(3);
	Teuchos::Array<double> ref_point(3);
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

	    TEST_ASSERT( cylinder.mapToReferenceFrame(point(),ref_point()) );
	    TEST_EQUALITY( ref_point[0], point[0] );
	    TEST_EQUALITY( ref_point[1], point[1] );
	    TEST_EQUALITY( ref_point[2], point[2] );

	    if ( centroid_distance <= radius + tol &&
		 centroid_z - length/2 <= point[2] + tol &&
		 centroid_z + length/2 >= point[2] - tol )
	    {
		TEST_ASSERT( cylinder.checkPointInclusion(tol,point()) );
	    }
	    else
	    {
		TEST_ASSERT( !cylinder.checkPointInclusion(tol,point()) );
	    }

	    cylinder.mapToPhysicalFrame(ref_point(),point());
	    TEST_EQUALITY( ref_point[0], point[0] );
	    TEST_EQUALITY( ref_point[1], point[1] );
	    TEST_EQUALITY( ref_point[2], point[2] );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstCylinder.cpp
//---------------------------------------------------------------------------//

