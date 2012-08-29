//---------------------------------------------------------------------------//
/*!
 * \file tstRendezvousMeshFromGeometry.cpp
 * \author Stuart R. Slattery
 * \brief RendezvousMesh unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_RendezvousMesh.hpp>
#include <DTK_GeometryTraits.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>

#include <MBInterface.hpp>
#include <MBRange.hpp>

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
TEUCHOS_UNIT_TEST( MeshContainer, box_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	double x_min = -(double) std::rand() / RAND_MAX;
	double y_min = -(double) std::rand() / RAND_MAX;
	double z_min = -(double) std::rand() / RAND_MAX;
	double x_max =  (double) std::rand() / RAND_MAX;
	double y_max =  (double) std::rand() / RAND_MAX;
	double z_max =  (double) std::rand() / RAND_MAX;
	boxes[i] = 
	    Box( x_min, y_min, z_min, x_max, y_max, z_max );
    }

    // Build a geometry manager.
    GeometryManager<Box> geometry_manager( boxes, comm, 3 );
    
    // Build a rendezvous mesh.
    Teuchos::RCP<RendezvousMesh<int> > mesh = 
	createRendezvousMeshFromGeometry<int,Box>( geometry_manager );

    // Get the moab interface.
    RendezvousMesh<int>::RCP_Moab moab = mesh->getMoab();
    
    // Check that it made the right number of hexahedrons.
    moab::Range mesh_elements;
    moab::ErrorCode error;
    error = moab()->get_entities_by_type( 0, moab::MBHEX, mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) mesh_elements.size() == num_boxes );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MeshContainer, cylinder_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
    }

    // Build a geometry manager.
    GeometryManager<Cylinder> geometry_manager( cylinders, comm, 3 );
    
    // Build a rendezvous mesh.
    Teuchos::RCP<RendezvousMesh<int> > mesh = 
	createRendezvousMeshFromGeometry<int,Cylinder>( geometry_manager );

    // Get the moab interface.
    RendezvousMesh<int>::RCP_Moab moab = mesh->getMoab();
    moab->write_mesh("out.vtk");    
    // Check that it made the right number of hexahedrons.
    moab::Range mesh_elements;
    moab::ErrorCode error;
    error = moab()->get_entities_by_type( 0, moab::MBHEX, mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) mesh_elements.size() == num_cylinders );
}

//---------------------------------------------------------------------------//
// end tstRendezvousMeshFromGeometry.cpp
//---------------------------------------------------------------------------//

