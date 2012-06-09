//---------------------------------------------------------------------------//
/*!
 * \file tstMeshContainer.cpp
 * \author Stuart R. Slattery
 * \brief MeshContainer unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MeshContainer.hpp>
#include <DTK_RendezvousMesh.hpp>
#include <DTK_CoreTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_TypeTraits.hpp>

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
// Mesh contianer creation funciton.
//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> buildMeshContainer()
{
    using namespace DataTransferKit;

    // Make some nodes.
    std::set<int> node_handles;
    std::vector<double> coords;

    node_handles.insert( 0 );
    coords.push_back( 0.0 ); coords.push_back( 0.0 ); coords.push_back( 0.0 );

    node_handles.insert( 1 );
    coords.push_back( 1.0 ); coords.push_back( 0.0 ); coords.push_back( 0.0 );

    node_handles.insert( 2 );
    coords.push_back( 1.0 ); coords.push_back( 1.0 ); coords.push_back( 0.0 );

    node_handles.insert( 3 );
    coords.push_back( 0.0 ); coords.push_back( 1.0 ); coords.push_back( 0.0 );

    node_handles.insert( 4 );
    coords.push_back( 0.0 ); coords.push_back( 0.0 ); coords.push_back( 1.0 );

    node_handles.insert( 5 );
    coords.push_back( 1.0 ); coords.push_back( 0.0 ); coords.push_back( 1.0 );

    node_handles.insert( 6 );
    coords.push_back( 1.0 ); coords.push_back( 1.0 ); coords.push_back( 1.0 );

    node_handles.insert( 7 );
    coords.push_back( 0.0 ); coords.push_back( 1.0 ); coords.push_back( 1.0 );

    node_handles.insert( 8 );
    coords.push_back( 0.0 ); coords.push_back( 0.0 ); coords.push_back( 2.0 );

    node_handles.insert( 9 );
    coords.push_back( 1.0 ); coords.push_back( 0.0 ); coords.push_back( 2.0 );

    node_handles.insert( 10 );
    coords.push_back( 1.0 ); coords.push_back( 1.0 ); coords.push_back( 2.0 );

    node_handles.insert( 11 );
    coords.push_back( 0.0 ); coords.push_back( 1.0 ); coords.push_back( 2.0 );

    // Make 2 hexahedrons.
    std::set<int> hex_handles;
    std::vector<int> hex_connectivity;
    
    hex_handles.insert( 0 );
    hex_connectivity.push_back( 0 ); hex_connectivity.push_back( 1 ); 
    hex_connectivity.push_back( 2 ); hex_connectivity.push_back( 3 ); 
    hex_connectivity.push_back( 4 ); hex_connectivity.push_back( 5 ); 
    hex_connectivity.push_back( 6 ); hex_connectivity.push_back( 7 ); 

    hex_handles.insert( 1 );
    hex_connectivity.push_back( 4 ); hex_connectivity.push_back( 5 ); 
    hex_connectivity.push_back( 6 ); hex_connectivity.push_back( 7 ); 
    hex_connectivity.push_back( 8 ); hex_connectivity.push_back( 9 ); 
    hex_connectivity.push_back( 10 ); hex_connectivity.push_back( 11 ); 

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );
    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );
    
    return MeshContainer<int>( node_handles, coords_array,
			       DTK_REGION, DTK_HEXAHEDRON, 8,
			       hex_handles, connectivity_array );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( MeshContainer, mesh_container_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshTraits< MeshContainer<int> > MT;
    MeshContainer<int> mesh_container = buildMeshContainer();
    Teuchos::RCP< RendezvousMesh<MT::handle_type> > mesh = 
	createRendezvousMesh( mesh_container );

    // Get the moab interface.
    moab::ErrorCode error;
    RendezvousMesh<MT::handle_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements = mesh->getElements();

    // Check the moab mesh element data.
    moab::Range::const_iterator element_iterator;
    MT::handle_type native_handle = 0;
    for ( element_iterator = mesh_elements.begin();
	  element_iterator != mesh_elements.end();
	  ++element_iterator, ++native_handle )
    {
	TEST_ASSERT( mesh->getNativeHandle( *element_iterator ) == 
		     native_handle );

	TEST_ASSERT( moab->type_from_handle( *element_iterator ) ==
		     moab::MBHEX );
    }

    // Check the moab mesh vertex data.
    moab::Range connectivity;
    error = moab->get_connectivity( mesh_elements, connectivity );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    std::vector<double> vertex_coords( 3 * connectivity.size() );
    error = moab->get_coords( connectivity, &vertex_coords[0] );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    std::vector<double>::const_iterator moab_coord_iterator = 
	vertex_coords.begin();
    typename MT::const_coordinate_iterator coord_iterator;

    for ( coord_iterator = MT::coordsBegin( mesh_container );
	  coord_iterator != MT::coordsEnd( mesh_container );
	  ++coord_iterator, ++moab_coord_iterator )
    {
	TEST_ASSERT( *coord_iterator == *moab_coord_iterator );
    }
}

//---------------------------------------------------------------------------//
// end tstMeshContainer.cpp
//---------------------------------------------------------------------------//
