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
#include <DTK_GeometricEntity.hpp>
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
// Tests
//---------------------------------------------------------------------------//
// Default constructor 1d test.
TEUCHOS_UNIT_TEST( Point, default_1d_constructor_test )
{
    using namespace DataTransferKit;

    // Make point.
    double x = 3.2;
    Point point(  0, 0, x );

    // Check GeometricEntity data.
    TEST_EQUALITY( point.entityType(), "DTK Point" );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 1 );
    TEST_EQUALITY( point.parametricDimension(), 0 );

    // Check the bounds.
    Teuchos::ArrayView<const double> point_coords;
    point.getCoordinates( point_coords );
    TEST_EQUALITY( point_coords[0], x );

    point.centroid( point_coords );
    TEST_EQUALITY( point_coords[0], x );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );
}

//---------------------------------------------------------------------------//
// Default constructor 2d test.
TEUCHOS_UNIT_TEST( Point, default_2d_constructor_test )
{
    using namespace DataTransferKit;

    // Make point.
    double x = 3.2;
    double y = -9.233;
    Point point(  0, 0, x, y );

    // Check GeometricEntity data.
    TEST_EQUALITY( point.entityType(), "DTK Point" );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 2 );
    TEST_EQUALITY( point.parametricDimension(), 0 );

    // Check the bounds.
    Teuchos::ArrayView<const double> point_coords;
    point.getCoordinates( point_coords );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );

    point.centroid( point_coords );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );
}

//---------------------------------------------------------------------------//
// Default constructor 3dtest.
TEUCHOS_UNIT_TEST( Point, default_3d_constructor_test )
{
    using namespace DataTransferKit;

    // Make point.
    double x = 3.2;
    double y = -9.233;
    double z = 1.3;
    Point point(  0, 0, x, y, z );

    // Check GeometricEntity data.
    TEST_EQUALITY( point.entityType(), "DTK Point" );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 3 );
    TEST_EQUALITY( point.parametricDimension(), 0 );

    // Check the bounds.
    Teuchos::ArrayView<const double> point_coords;
    point.getCoordinates( point_coords );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );

    point.centroid( point_coords );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );
}

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

    // Check GeometricEntity data.
    TEST_EQUALITY( point.entityType(), "DTK Point" );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 1 );
    TEST_EQUALITY( point.parametricDimension(), 0 );

    // Check the bounds.
    Teuchos::ArrayView<const double> point_coords;
    point.getCoordinates( point_coords );
    TEST_EQUALITY( point_coords[0], x );

    point.centroid( point_coords );
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

    // Check GeometricEntity data.
    TEST_EQUALITY( point.entityType(), "DTK Point" );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 2 );
    TEST_EQUALITY( point.parametricDimension(), 0 );

    // Check the bounds.
    Teuchos::ArrayView<const double> point_coords;
    point.getCoordinates( point_coords );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );

    point.centroid( point_coords );
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

    // Check GeometricEntity data.
    TEST_EQUALITY( point.entityType(), "DTK Point" );
    TEST_EQUALITY( point.id(), 0 );
    TEST_EQUALITY( point.ownerRank(), 0 );
    TEST_EQUALITY( point.physicalDimension(), 3 );
    TEST_EQUALITY( point.parametricDimension(), 0 );

    // Check the bounds.
    Teuchos::ArrayView<const double> point_coords;
    point.getCoordinates( point_coords );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );

    point.centroid( point_coords );
    TEST_EQUALITY( point_coords[0], x );
    TEST_EQUALITY( point_coords[1], y );
    TEST_EQUALITY( point_coords[2], z );

    // Compute the measure.
    TEST_EQUALITY( point.measure(), 0.0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Point, communication_test )
{
    using namespace DataTransferKit;

    // Register the point class with the geometric entity class so that we can
    // send the point to a different rank and reconstruct it.
    GeometricEntity::setDerivedClassFactory(
	Teuchos::abstractFactoryStd<GeometricEntity,Point>() );

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Make a point.
    double x = 3.2;
    double y = -9.233;
    double z = 1.3;

    Teuchos::RCP<GeometricEntity> entity;
    if ( 0 == comm_rank )
    {
	entity = Teuchos::rcp( new Point(0, 0, x, y, z) );
    }

    // Broadcast the point with indirect serialization through the geometric
    // entity api.
    Teuchos::broadcast( *comm, 0, 
			Teuchos::Ptr<Teuchos::RCP<GeometricEntity> >(&entity) );

    // Check the coordinates.
    Teuchos::ArrayView<const double> coords;
    entity->centroid( coords );
    TEST_ASSERT( coords[0] == x );
    TEST_ASSERT( coords[1] == y );
    TEST_ASSERT( coords[2] == z );

    // Broadcast an array.
    Teuchos::Array<Teuchos::RCP<GeometricEntity> > points(2);
    double x_1 = 3.2 + comm_size;
    double y_1 = -9.233 + comm_size;
    double z_1 = 1.3 + comm_size;
    double x_2 = 3.2 - comm_size;
    double y_2 = -9.233 - comm_size;
    double z_2 = 1.3 - comm_size;
    if ( 0 == comm->getRank() )
    {
	points[0] = Teuchos::rcp( new Point(0, 0, x_1, y_1, z_1) );
	points[1] = Teuchos::rcp( new Point(1, 0, x_2, y_2, z_2) );
    }
    Teuchos::broadcast( *comm, 0, points() );
    points[0]->centroid( coords );
    TEST_EQUALITY( coords[0], x_1 );
    TEST_EQUALITY( coords[1], y_1 );
    TEST_EQUALITY( coords[2], z_1 );
    points[1]->centroid( coords );
    TEST_EQUALITY( coords[0], x_2 );
    TEST_EQUALITY( coords[1], y_2 );
    TEST_EQUALITY( coords[2], z_2 );
}

//---------------------------------------------------------------------------//
// end tstPoint.cpp
//---------------------------------------------------------------------------//

