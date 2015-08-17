//---------------------------------------------------------------------------//
/*!
 * \file tstMeshManager.cpp
 * \author Stuart R. Slattery
 * \brief Mesh manager unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
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
// Mesh container creation functions.
//---------------------------------------------------------------------------//
// Line segment mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildLineContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> line_handles;
    Teuchos::Array<unsigned long int> line_connectivity;
    
    // handles
    line_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	line_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> line_handle_array( line_handles.size() );
    std::copy( line_handles.begin(), line_handles.end(), 
	       line_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( line_connectivity.size() );
    std::copy( line_connectivity.begin(), line_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_LINE_SEGMENT, num_vertices,
					      line_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tri mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildTriContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> tri_handles;
    Teuchos::Array<unsigned long int> tri_connectivity;
    
    // handles
    tri_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	tri_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> tri_handle_array( tri_handles.size() );
    std::copy( tri_handles.begin(), tri_handles.end(), 
	       tri_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( tri_connectivity.size() );
    std::copy( tri_connectivity.begin(), tri_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_TRIANGLE, num_vertices,
					      tri_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Quad mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildQuadContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> quad_handles;
    Teuchos::Array<unsigned long int> quad_connectivity;
    
    // handles
    quad_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	quad_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_QUADRILATERAL, num_vertices,
					      quad_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tet mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildTetContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> tet_handles;
    Teuchos::Array<unsigned long int> tet_connectivity;
    
    // handles
    tet_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	tet_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> tet_handle_array( tet_handles.size() );
    std::copy( tet_handles.begin(), tet_handles.end(), 
	       tet_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( tet_connectivity.size() );
    std::copy( tet_connectivity.begin(), tet_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_TETRAHEDRON, num_vertices,
					      tet_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Hex mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildHexContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> hex_handles;
    Teuchos::Array<unsigned long int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_HEXAHEDRON, num_vertices,
					      hex_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildParallelHexContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> hex_handles;
    Teuchos::Array<unsigned long int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_HEXAHEDRON, num_vertices,
					      hex_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildPyramidContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> pyramid_handles;
    Teuchos::Array<unsigned long int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_PYRAMID, num_vertices,
					      pyramid_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<unsigned long int> > buildWedgeContainer()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<unsigned long int> vertex_handles;
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
    Teuchos::Array<unsigned long int> wedge_handles;
    Teuchos::Array<unsigned long int> wedge_connectivity;
    
    // handles
    wedge_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	wedge_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> wedge_handle_array( wedge_handles.size() );
    std::copy( wedge_handles.begin(), wedge_handles.end(), 
	       wedge_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array( wedge_connectivity.size() );
    std::copy( wedge_connectivity.begin(), wedge_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<unsigned long int>( vertex_dim, vertex_handle_array, coords_array,
					      DTK_WEDGE, num_vertices,
					      wedge_handle_array, connectivity_array,
					      permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshContainer, line_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildLineContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 1 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 1 );

    // Mesh parameters.
    int num_vertices = 2;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Basic container info.
	TEST_ASSERT( Tools::numElements( *(*block_iterator) ) == 1 );
	TEST_ASSERT( (int) Tools::numVertices( *(*block_iterator) ) == num_vertices );

	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tri_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTriContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Mesh parameters.
    int num_vertices = 3;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	TEST_ASSERT( Tools::numElements( *(*block_iterator) ) == 1 );
	TEST_ASSERT( (int) Tools::numVertices( *(*block_iterator) ) == num_vertices );

	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 

	// y
	TEST_ASSERT( coords_view[3] == 0.0 ); 
	TEST_ASSERT( coords_view[4] == 0.0 ); 
	TEST_ASSERT( coords_view[5] == 1.0 ); 

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshContainer, quad_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildQuadContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Mesh parameters.
    int num_vertices = 4;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );

	// y
	TEST_ASSERT( coords_view[4] == 0.0 ); 
	TEST_ASSERT( coords_view[5] == 0.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 1.0 ); 

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding boxes.
    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Tet mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tet_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTetContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 4;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );

	// y
	TEST_ASSERT( coords_view[4] == 0.0 ); 
	TEST_ASSERT( coords_view[5] == 0.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 1.0 ); 

	// z
	TEST_ASSERT( coords_view[8]  == 0.0 );
	TEST_ASSERT( coords_view[9]  == 0.0 );
	TEST_ASSERT( coords_view[10] == 0.0 );
	TEST_ASSERT( coords_view[11] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, hex_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildHexContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 8;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 0.0 );
	TEST_ASSERT( coords_view[5] == 1.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 0.0 ); 

	// y
	TEST_ASSERT( coords_view[8]  == 0.0 ); 
	TEST_ASSERT( coords_view[9]  == 0.0 ); 
	TEST_ASSERT( coords_view[10] == 1.0 ); 
	TEST_ASSERT( coords_view[11] == 1.0 ); 
	TEST_ASSERT( coords_view[12] == 0.0 ); 
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 1.0 );
	TEST_ASSERT( coords_view[15] == 1.0 );

	// z
	TEST_ASSERT( coords_view[16] == 0.0 );
	TEST_ASSERT( coords_view[17] == 0.0 );
	TEST_ASSERT( coords_view[18] == 0.0 );
	TEST_ASSERT( coords_view[19] == 0.0 );
	TEST_ASSERT( coords_view[20] == 1.0 );
	TEST_ASSERT( coords_view[21] == 1.0 );
	TEST_ASSERT( coords_view[22] == 1.0 );
	TEST_ASSERT( coords_view[23] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
TEUCHOS_UNIT_TEST( MeshContainer, pyramid_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildPyramidContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 5;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 0.5 );

	// y
	TEST_ASSERT( coords_view[5]  == 0.0 ); 
	TEST_ASSERT( coords_view[6]  == 0.0 ); 
	TEST_ASSERT( coords_view[7] == 1.0 ); 
	TEST_ASSERT( coords_view[8] == 1.0 ); 
	TEST_ASSERT( coords_view[9] == 0.5 ); 

	// z
	TEST_ASSERT( coords_view[10] == 0.0 );
	TEST_ASSERT( coords_view[11] == 0.0 );
	TEST_ASSERT( coords_view[12] == 0.0 );
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
TEUCHOS_UNIT_TEST( MeshContainer, wedge_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildWedgeContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 6;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 0.5 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 1.0 );
	TEST_ASSERT( coords_view[5] == 0.5 );

	// y
	TEST_ASSERT( coords_view[6]  == 0.0 ); 
	TEST_ASSERT( coords_view[7]  == 0.0 ); 
	TEST_ASSERT( coords_view[8]  == 1.0 ); 
	TEST_ASSERT( coords_view[9]  == 0.0 ); 
	TEST_ASSERT( coords_view[10] == 0.0 ); 
	TEST_ASSERT( coords_view[11] == 1.0 ); 

	// z
	TEST_ASSERT( coords_view[12] == 0.0 );
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 0.0 );
	TEST_ASSERT( coords_view[15] == 1.0 );
	TEST_ASSERT( coords_view[16] == 1.0 );
	TEST_ASSERT( coords_view[17] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, parallel_hex_manager_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildParallelHexContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 8;

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    Tools::coordsView( *(*block_iterator) );
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 0.0 );
	TEST_ASSERT( coords_view[5] == 1.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 0.0 ); 

	// y
	TEST_ASSERT( coords_view[8]  == 0.0 ); 
	TEST_ASSERT( coords_view[9]  == 0.0 ); 
	TEST_ASSERT( coords_view[10] == 1.0 ); 
	TEST_ASSERT( coords_view[11] == 1.0 ); 
	TEST_ASSERT( coords_view[12] == 0.0 ); 
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 1.0 );
	TEST_ASSERT( coords_view[15] == 1.0 );

	// z
	TEST_ASSERT( coords_view[16] == my_rank );
	TEST_ASSERT( coords_view[17] == my_rank );
	TEST_ASSERT( coords_view[18] == my_rank );
	TEST_ASSERT( coords_view[19] == my_rank );
	TEST_ASSERT( coords_view[20] == my_rank+1 );
	TEST_ASSERT( coords_view[21] == my_rank+1 );
	TEST_ASSERT( coords_view[22] == my_rank+1 );
	TEST_ASSERT( coords_view[23] == my_rank+1 );

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == my_size );
}

//---------------------------------------------------------------------------//
// 2d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 2d_hybrid_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    mesh_blocks[0] = buildTriContainer();
    mesh_blocks[1] = buildQuadContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 2 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );

    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	TEST_ASSERT( Tools::numElements( *(*block_iterator) ) == 1 );
	int num_vertices = (int) Tools::numVertices( *(*block_iterator) );

	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// 3d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 3d_hybrid_manager_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<unsigned long int> MeshType;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 3 );
    mesh_blocks[0] = buildTetContainer();
    mesh_blocks[1] = buildHexContainer();
    mesh_blocks[2] = buildPyramidContainer();

    // Create a mesh manager.
    MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 3 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Check the mesh data.
    MeshManager<MeshType>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	TEST_ASSERT( Tools::numElements( *(*block_iterator) ) == 1 );
	int num_vertices = (int) Tools::numVertices( *(*block_iterator) );

	// Vertices.
	Teuchos::ArrayRCP<const unsigned long int> vertices_view = 
	    Tools::verticesView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Elements.
	Teuchos::ArrayRCP<const unsigned long int> elements_view =
	    Tools::elementsView( *(*block_iterator) );
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
	    Tools::connectivityView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    Tools::permutationView( *(*block_iterator) );
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// end tstMeshManager.cpp
//---------------------------------------------------------------------------//

