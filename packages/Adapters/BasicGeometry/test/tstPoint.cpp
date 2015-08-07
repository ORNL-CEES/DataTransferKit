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
 * \file tstPoint.cpp
 * \author Stuart R. Slattery
 * \brief Point unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_Point.hpp>
#include <DTK_Entity.hpp>
#include <DTK_Box.hpp>

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
// Tests
//---------------------------------------------------------------------------//
// Array constructor 1d test.
TEUCHOS_UNIT_TEST( Point, array_1d_constructor_test )
{
    using namespace DataTransferKit;

    // Make point.
    double x = 3.2;
    Teuchos::Array<double> p(1);
    p[0] = x;
    Point point(  0, 0, p );

    // Check Entity data.
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.topologicalDimension(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 1 );

    // Check the bounds.
    Teuchos::Array<double> point_coords(1);
    point.getCoordinates( point_coords() );
    TEST_EQUALITY( point_coords[0], x );

    point.centroid( point_coords() );
    TEST_EQUALITY( point_coords[0], x );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );
}

//---------------------------------------------------------------------------//
// Array constructor 2d test.
TEUCHOS_UNIT_TEST( Point, array_2d_constructor_test )
{
    using namespace DataTransferKit;

    // Make point.
    double x = 3.2;
    double y = -9.233;
    Teuchos::Array<double> p(2);
    p[0] = x;
    p[1] = y;
    Point point(  0, 0, p );

    // Check Entity data.
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.topologicalDimension(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 2 );

    // Check the bounds.
    Teuchos::Array<double> point_coords(2);
    point.getCoordinates( point_coords() );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );

    point.centroid( point_coords() );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );
}

//---------------------------------------------------------------------------//
// Array constructor 3d test.
TEUCHOS_UNIT_TEST( Point, array_3d_constructor_test )
{
    using namespace DataTransferKit;

    // Make point.
    double x = 3.2;
    double y = -9.233;
    double z = 1.3;
    Teuchos::Array<double> p(3);
    p[0] = x;
    p[1] = y;
    p[2] = z;
    Point point(  0, 0, p );

    // Check Entity data.
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.topologicalDimension(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 3 );

    // Check the coordinates.
    Teuchos::Array<double> point_coords(3);
    point.getCoordinates( point_coords() );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );
    point_coords.clear();
    point_coords.resize(3);
    point.centroid( point_coords() );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );

    // Copy to the BasicGeometryEntity interface and check.
    BasicGeometryEntity geom_entity = point;
    TEST_EQUALITY( geom_entity.id(), 0 );
    TEST_EQUALITY( geom_entity.ownerRank(), 0 );
    TEST_EQUALITY( geom_entity.topologicalDimension(), 0 );
    TEST_EQUALITY( geom_entity.physicalDimension(), 3 );
    TEST_EQUALITY( geom_entity.measure(), 0.0 );
    point_coords.clear();
    point_coords.resize(3);
    geom_entity.centroid( point_coords() );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );
    Teuchos::Array<double> ref_point(3);
    geom_entity.mapToReferenceFrame( p(), ref_point() );
    TEST_EQUALITY( ref_point[0], x );
    TEST_EQUALITY( ref_point[1], y );
    TEST_EQUALITY( ref_point[2], z );
    TEST_ASSERT( geom_entity.checkPointInclusion(1.0e-6,p()) );
    p[0] += 1.0e-7;
    TEST_ASSERT( geom_entity.checkPointInclusion(1.0e-6,p()) );
    p[0] += 1.0e-5;
    TEST_ASSERT( !geom_entity.checkPointInclusion(1.0e-6,p()) );
    geom_entity.mapToPhysicalFrame( ref_point(), p() );
    TEST_EQUALITY( p[0], x );
    TEST_EQUALITY( p[1], y );
    TEST_EQUALITY( p[2], z );

    // Copy to the Entity interface and check.
    Entity entity = point;
    TEST_EQUALITY( entity.id(), 0 );
    TEST_EQUALITY( entity.ownerRank(), 0 );
    TEST_EQUALITY( geom_entity.topologicalDimension(), 0 );
    TEST_EQUALITY( entity.physicalDimension(), 3 );
}

//---------------------------------------------------------------------------//
// end tstPoint.cpp
//---------------------------------------------------------------------------//

