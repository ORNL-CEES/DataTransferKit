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
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
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
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    // handles
    for ( int i = 0; i < 12; ++i )
    {
	node_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 2.0 );
    coords.push_back( 2.0 );
    coords.push_back( 2.0 );
    coords.push_back( 2.0 );

    // Make 2 hexahedrons.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 0 );
    hex_handles.push_back( 1 );

    // 0
    hex_connectivity.push_back( 0 );
    hex_connectivity.push_back( 4 ); 

    // 1
    hex_connectivity.push_back( 1 ); 
    hex_connectivity.push_back( 5 );  

    // 2
    hex_connectivity.push_back( 2 );
    hex_connectivity.push_back( 6 ); 

    // 3
    hex_connectivity.push_back( 3 ); 
    hex_connectivity.push_back( 7 ); 

    // 4
    hex_connectivity.push_back( 4 );
    hex_connectivity.push_back( 8 ); 
   
    // 5
    hex_connectivity.push_back( 5 ); 
    hex_connectivity.push_back( 9 ); 

    // 6
    hex_connectivity.push_back( 6 ); 
    hex_connectivity.push_back( 10 ); 

    // 7
    hex_connectivity.push_back( 7 ); 
    hex_connectivity.push_back( 11 ); 

    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );
    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );
    Teuchos::ArrayRCP<int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );
    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );
    Teuchos::ArrayRCP<std::size_t> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( 3, node_handle_array, coords_array,
			       DTK_HEXAHEDRON, 8,
			       hex_handle_array, connectivity_array,
			       permutation_list );
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
    Teuchos::RCP< RendezvousMesh<MT::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_container );

    // Get the moab interface.
    moab::ErrorCode error;
    RendezvousMesh<MT::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements = mesh->getElements();

    // Check the moab mesh element data.
    moab::Range::const_iterator element_iterator;
    MT::global_ordinal_type native_handle = 0;
    for ( element_iterator = mesh_elements.begin();
	  element_iterator != mesh_elements.end();
	  ++element_iterator, ++native_handle )
    {
	TEST_ASSERT( mesh->getNativeOrdinal( *element_iterator ) == 
		     native_handle );

	TEST_ASSERT( moab->type_from_handle( *element_iterator ) ==
		     moab::MBHEX );
    }

    // Check the moab mesh vertex data.
    moab::Range connectivity;
    error = moab->get_connectivity( mesh_elements, connectivity );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    Teuchos::Array<double> vertex_coords( 3 * connectivity.size() );
    error = moab->get_coords( connectivity, &vertex_coords[0] );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    int num_nodes = connectivity.size();
    Teuchos::Array<double>::const_iterator moab_coord_iterator;
    typename MT::const_coordinate_iterator coord_iterator = 
	MT::coordsBegin( mesh_container );
    int i = 0;
    for ( moab_coord_iterator = vertex_coords.begin();
	  moab_coord_iterator != vertex_coords.end(); ++i )
    {
	for ( int d = 0; d < 3; ++d, ++moab_coord_iterator )
	{
	    TEST_ASSERT( coord_iterator[d*num_nodes + i] == 
			 *moab_coord_iterator );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstMeshContainer.cpp
//---------------------------------------------------------------------------//
