//---------------------------------------------------------------------------//
/*!
 * \file tstGeometryRendezvous.cpp
 * \author Stuart R. Slattery
 * \brief Unit tests for recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_Partitioner.hpp>
#include <DTK_PartitionerFactory.hpp>
#include <DTK_GeometryRendezvous.hpp>
#include <DTK_BoundingBox.hpp>
#include <DTK_GeometryTraits.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Global test variables.
//---------------------------------------------------------------------------//
int rand_size = 1000;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// cylinder test
TEUCHOS_UNIT_TEST( GeometryRendezvous, cylinder_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    Teuchos::ArrayRCP<int> gids( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
	gids[i] = i + my_rank*num_cylinders;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Cylinder,int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Cylinder,int>( 
			  cylinders, gids, comm, 3 ) );

    // Build a rendezvous.
    BoundingBox global_box( -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax() );
    GeometryRendezvous<Cylinder,int> rendezvous( 
	getDefaultComm<int>(), 3, global_box );
    rendezvous.build( geometry_manager );
    
    // Check the rendezvous point search.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
}

//---------------------------------------------------------------------------//
// box test
TEUCHOS_UNIT_TEST( GeometryRCB, box_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    Teuchos::ArrayRCP<int> gids( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	double x_min = -(double) std::rand() / RAND_MAX + my_rank;
	double y_min = -(double) std::rand() / RAND_MAX + my_rank;
	double z_min = -(double) std::rand() / RAND_MAX + my_rank;
	double x_max =  (double) std::rand() / RAND_MAX + my_rank;
	double y_max =  (double) std::rand() / RAND_MAX + my_rank;
	double z_max =  (double) std::rand() / RAND_MAX + my_rank;
	boxes[i] = 
	    Box( x_min, y_min, z_min, x_max, y_max, z_max );
	gids[i] = i + my_rank*num_boxes;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Box,int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Box,int>( 
			  boxes, gids, comm, 3 ) );

    // Build a rendezvous.
    BoundingBox global_box( -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax() );
    GeometryRendezvous<Box,int> rendezvous( 
	getDefaultComm<int>(), 3, global_box );
    rendezvous.build( geometry_manager );
    
    // Check the rendezvous point search.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
}

//---------------------------------------------------------------------------//
// cylinder test
TEUCHOS_UNIT_TEST( GeometryRendezvous, cylinder_unsigned_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random cylinders.
    unsigned int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    Teuchos::ArrayRCP<unsigned int> gids( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
	gids[i] = i + my_rank*num_cylinders;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Cylinder,unsigned int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Cylinder,unsigned int>( 
			  cylinders, gids, comm, 3 ) );

    // Build a rendezvous.
    BoundingBox global_box( -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax() );
    GeometryRendezvous<Cylinder,unsigned int> rendezvous( 
	getDefaultComm<int>(), 3, global_box );
    rendezvous.build( geometry_manager );
    
    // Check the rendezvous point search.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
}

//---------------------------------------------------------------------------//
// box test
TEUCHOS_UNIT_TEST( GeometryRCB, box_unsigned_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    Teuchos::ArrayRCP<unsigned int> gids( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	double x_min = -(double) std::rand() / RAND_MAX + my_rank;
	double y_min = -(double) std::rand() / RAND_MAX + my_rank;
	double z_min = -(double) std::rand() / RAND_MAX + my_rank;
	double x_max =  (double) std::rand() / RAND_MAX + my_rank;
	double y_max =  (double) std::rand() / RAND_MAX + my_rank;
	double z_max =  (double) std::rand() / RAND_MAX + my_rank;
	boxes[i] = 
	    Box( x_min, y_min, z_min, x_max, y_max, z_max );
	gids[i] = i + my_rank*num_boxes;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Box,unsigned int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Box,unsigned int>( 
			  boxes, gids, comm, 3 ) );

    // Build a rendezvous.
    BoundingBox global_box( -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax() );
    GeometryRendezvous<Box,unsigned int> rendezvous( 
	getDefaultComm<int>(), 3, global_box );
    rendezvous.build( geometry_manager );
    
    // Check the rendezvous point search.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
}

//---------------------------------------------------------------------------//
// end tstGeometryRendezvous.cpp
//---------------------------------------------------------------------------//
