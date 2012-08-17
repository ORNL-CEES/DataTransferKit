//---------------------------------------------------------------------------//
/*!
 * \file tstRendezvousMesh.cpp
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
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
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
// Mesh container creation functions.
//---------------------------------------------------------------------------//
// Line segment mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildLineContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 1;
    int num_vertices = 2;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 

    // Make the line.
    Teuchos::Array<int> line_handles;
    Teuchos::Array<int> line_connectivity;
    
    // handles
    line_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	line_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> line_handle_array( line_handles.size() );
    std::copy( line_handles.begin(), line_handles.end(), 
	       line_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( line_connectivity.size() );
    std::copy( line_connectivity.begin(), line_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_LINE_SEGMENT, num_vertices,
				line_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tri mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildTriContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 2;
    int num_vertices = 3;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 

    // Make the triahedron.
    Teuchos::Array<int> tri_handles;
    Teuchos::Array<int> tri_connectivity;
    
    // handles
    tri_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	tri_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> tri_handle_array( tri_handles.size() );
    std::copy( tri_handles.begin(), tri_handles.end(), 
	       tri_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( tri_connectivity.size() );
    std::copy( tri_connectivity.begin(), tri_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_TRIANGLE, num_vertices,
				tri_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Quad mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildQuadContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 2;
    int num_vertices = 4;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 

    // Make the quadahedron.
    Teuchos::Array<int> quad_handles;
    Teuchos::Array<int> quad_connectivity;
    
    // handles
    quad_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	quad_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_QUADRILATERAL, num_vertices,
				quad_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tet mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildTetContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 4;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );

    // Make the tetahedron.
    Teuchos::Array<int> tet_handles;
    Teuchos::Array<int> tet_connectivity;
    
    // handles
    tet_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	tet_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> tet_handle_array( tet_handles.size() );
    std::copy( tet_handles.begin(), tet_handles.end(), 
	       tet_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( tet_connectivity.size() );
    std::copy( tet_connectivity.begin(), tet_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_TETRAHEDRON, num_vertices,
				tet_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Hex mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildHexContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 8;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
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

    // y
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

    // Make the hexahedron.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_HEXAHEDRON, num_vertices,
				hex_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildParallelHexContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 8;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
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

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );

    // z
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );

    // Make the hexahedron.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_HEXAHEDRON, num_vertices,
				hex_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildPyramidContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 5;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 0.5 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.5 ); 

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );

    // Make the pyramidahedron.
    Teuchos::Array<int> pyramid_handles;
    Teuchos::Array<int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_PYRAMID, num_vertices,
				pyramid_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildWedgeContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 6;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.5 ); 
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 0.5 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 ); 

    // Make the wedge.
    Teuchos::Array<int> wedge_handles;
    Teuchos::Array<int> wedge_connectivity;
    
    // handles
    wedge_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	wedge_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> wedge_handle_array( wedge_handles.size() );
    std::copy( wedge_handles.begin(), wedge_handles.end(), 
	       wedge_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( wedge_connectivity.size() );
    std::copy( wedge_connectivity.begin(), wedge_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_WEDGE, num_vertices,
				wedge_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshContainer, line_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildLineContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 1 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 1 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = 
	mesh->getMoab();
    
    // Elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_LINE_SEGMENT ] ==
		     element_type );
	TEST_ASSERT( mesh->getNativeOrdinal( mesh_elements[i] ) == 12 );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == 
		 MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == 
			 mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tri_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTriContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_TRIANGLE ] ==
		     element_type );
	TEST_ASSERT( mesh->getNativeOrdinal( mesh_elements[i] ) == 12 );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == 
		 MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshContainer, quad_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildQuadContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = 
	mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_QUADRILATERAL ] ==
		     element_type );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// Tet mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tet_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTetContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_TETRAHEDRON ] ==
		     element_type );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == 
		 MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// Hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, hex_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildHexContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_HEXAHEDRON ] ==
		     element_type );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == 
		 MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
TEUCHOS_UNIT_TEST( MeshContainer, pyramid_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildPyramidContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_PYRAMID ] ==
		     element_type );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == 
		 MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// Wedge mesh.
TEUCHOS_UNIT_TEST( MeshContainer, wedge_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildWedgeContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_WEDGE ] ==
		     element_type );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == 
		 MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, parallel_hex_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildParallelHexContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );
    for ( int i = 0; i < (int) mesh_elements.size(); ++i )
    {
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[i] );
	TEST_ASSERT( moab_topology_table[ DTK_HEXAHEDRON ] ==
		     element_type );
    }

    // Vertices
    moab::Range vertices;
    error = moab->get_connectivity( mesh_elements, vertices );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    TEST_ASSERT( (int) vertices.size() == 
		 MT::verticesPerElement( *mesh_blocks[0] ) );

    // Coords.
    int vertex_dim = MT::vertexDim( *mesh_blocks[0] );
    std::vector<double> mb_coords( 3*vertices.size() );
    error = moab->get_coords( vertices, &mb_coords[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    Teuchos::ArrayRCP<const double> coords_view = 
	Tools::coordsView( *mesh_blocks[0] );
    for ( int i = 0; i < (int) vertices.size(); ++i )
    {
	for ( int d = 0; d < vertex_dim; ++d )
	{
	    TEST_ASSERT( coords_view[vertices.size()*d + i] == mb_coords[3*i+d] ); 
	}
    }
}

//---------------------------------------------------------------------------//
// 2d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 2d_hybrid_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    mesh_blocks[0] = buildTriContainer();
    mesh_blocks[1] = buildQuadContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 2 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    moab::EntityHandle root_set = moab->get_root_set();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	int block_id = std::distance( mesh_manager.blocksBegin(),
				      block_iterator );
	// Elements.
	int block_topology = MT::elementTopology( *(*block_iterator) );
	Teuchos::ArrayRCP<const int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[block_id] );
	TEST_ASSERT( moab_topology_table[ block_topology ] ==
		     element_type );

	// Vertices.
	int num_vertices = MT::verticesPerElement( *(*block_iterator) );
	moab::Range block_elements;
	error = moab->get_entities_by_type( 
	    root_set, element_type, block_elements );
	TEST_ASSERT( error == moab::MB_SUCCESS );
	moab::Range vertices;
	error = moab->get_connectivity( block_elements, vertices );
	TEST_ASSERT( error == moab::MB_SUCCESS );
	TEST_ASSERT( (int) vertices.size() == num_vertices );

	// Coords.
	int vertex_dim = MT::vertexDim( *(*block_iterator) );
	std::vector<double> mb_coords( 3*vertices.size() );
	error = moab->get_coords( vertices, &mb_coords[0] );
	TEST_ASSERT( error == moab::MB_SUCCESS );

	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	for ( int i = 0; i < (int) vertices.size(); ++i )
	{
	    for ( int d = 0; d < vertex_dim; ++d )
	    {
		TEST_ASSERT( coords_view[vertices.size()*d + i] == 
			     mb_coords[3*i+d] ); 
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// 3d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 3d_hybrid_rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 3 );
    mesh_blocks[0] = buildTetContainer();
    mesh_blocks[1] = buildPyramidContainer();
    mesh_blocks[2] = buildHexContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 3 );
    TEST_ASSERT( mesh_manager.comm() == getDefaultComm<int>() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Create a rendezvous mesh.
    moab::ErrorCode error;
    Teuchos::RCP< RendezvousMesh<MeshType::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MeshType::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    moab::EntityHandle root_set = moab->get_root_set();
    
    // Grab the elements.
    moab::Range mesh_elements;
    error = moab()->get_entities_by_dimension( 0, mesh_manager.dim(), 
					       mesh_elements );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    TEST_ASSERT( (int) mesh_elements.size() == mesh_manager.localNumElements() );

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	int block_id = std::distance( mesh_manager.blocksBegin(),
				      block_iterator );

	// Elements.
	int block_topology = MT::elementTopology( *(*block_iterator) );
	Teuchos::ArrayRCP<const int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );
	moab::EntityType element_type = 
	    moab->type_from_handle( mesh_elements[block_id] );
	TEST_ASSERT( moab_topology_table[ block_topology ] ==
		     element_type );

	// Vertices.
	int num_vertices = MT::verticesPerElement( *(*block_iterator) );
	moab::Range block_elements;
	error = moab->get_entities_by_type( 
	    root_set, element_type, block_elements );
	TEST_ASSERT( error == moab::MB_SUCCESS );
	moab::Range vertices;
	error = moab->get_connectivity( block_elements, vertices );
	TEST_ASSERT( error == moab::MB_SUCCESS );
	TEST_ASSERT( (int) vertices.size() == num_vertices );

	// Coords.
	int vertex_dim = MT::vertexDim( *(*block_iterator) );
	std::vector<double> mb_coords( 3*vertices.size() );
	error = moab->get_coords( vertices, &mb_coords[0] );
	TEST_ASSERT( error == moab::MB_SUCCESS );

	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	for ( int i = 0; i < (int) vertices.size(); ++i )
	{
	    for ( int d = 0; d < vertex_dim; ++d )
	    {
		TEST_ASSERT( coords_view[vertices.size()*d + i] == 
			     mb_coords[3*i+d] ); 
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// end tstRendezvousMesh.cpp
//---------------------------------------------------------------------------//

