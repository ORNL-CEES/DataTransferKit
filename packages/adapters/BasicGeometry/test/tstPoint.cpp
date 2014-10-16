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
    TEST_EQUALITY( point.entityType(), ENTITY_TYPE_NODE );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
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
    TEST_EQUALITY( point.entityType(), ENTITY_TYPE_NODE );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
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
    TEST_EQUALITY( point.entityType(), ENTITY_TYPE_NODE );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 3 );

    // Check the bounds.
    Teuchos::Array<double> point_coords(3);
    point.getCoordinates( point_coords() );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );

    point.centroid( point_coords() );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );
}

//---------------------------------------------------------------------------//
// end tstPoint.cpp
//---------------------------------------------------------------------------//

