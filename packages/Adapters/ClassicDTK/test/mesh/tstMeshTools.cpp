//---------------------------------------------------------------------------//
/*!
 * \file tstMeshTools.cpp
 * \author Stuart R. Slattery
 * \brief Mesh tools unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_BoundingBox.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTypes.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template <class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Mesh container creation functions.
//---------------------------------------------------------------------------//
// Line segment mesh.
DataTransferKit::MeshContainer<unsigned long int> buildLineContainer()
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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> line_handle_array(
        line_handles.size() );
    std::copy( line_handles.begin(), line_handles.end(),
               line_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        line_connectivity.size() );
    std::copy( line_connectivity.begin(), line_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_LINE_SEGMENT,
        num_vertices, line_handle_array, connectivity_array, permutation_list );
}

//---------------------------------------------------------------------------//
// Tri mesh.
DataTransferKit::MeshContainer<unsigned long int> buildTriContainer()
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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> tri_handle_array( tri_handles.size() );
    std::copy( tri_handles.begin(), tri_handles.end(),
               tri_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        tri_connectivity.size() );
    std::copy( tri_connectivity.begin(), tri_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_TRIANGLE,
        num_vertices, tri_handle_array, connectivity_array, permutation_list );
}

//---------------------------------------------------------------------------//
// Quad mesh.
DataTransferKit::MeshContainer<unsigned long int> buildQuadContainer()
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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> quad_handle_array(
        quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(),
               quad_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_QUADRILATERAL,
        num_vertices, quad_handle_array, connectivity_array, permutation_list );
}

//---------------------------------------------------------------------------//
// Tet mesh.
DataTransferKit::MeshContainer<unsigned long int> buildTetContainer()
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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> tet_handle_array( tet_handles.size() );
    std::copy( tet_handles.begin(), tet_handles.end(),
               tet_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        tet_connectivity.size() );
    std::copy( tet_connectivity.begin(), tet_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_TETRAHEDRON,
        num_vertices, tet_handle_array, connectivity_array, permutation_list );
}

//---------------------------------------------------------------------------//
// Hex mesh.
DataTransferKit::MeshContainer<unsigned long int> buildHexContainer()
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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(),
               hex_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_HEXAHEDRON,
        num_vertices, hex_handle_array, connectivity_array, permutation_list );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
DataTransferKit::MeshContainer<unsigned long int> buildWedgeContainer()
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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> wedge_handle_array(
        wedge_handles.size() );
    std::copy( wedge_handles.begin(), wedge_handles.end(),
               wedge_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        wedge_connectivity.size() );
    std::copy( wedge_connectivity.begin(), wedge_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_WEDGE, num_vertices,
        wedge_handle_array, connectivity_array, permutation_list );
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
DataTransferKit::MeshContainer<unsigned long int> buildParallelHexContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<unsigned long int>()->getRank();

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
    coords.push_back( my_rank + 1 );
    coords.push_back( my_rank + 1 );
    coords.push_back( my_rank + 1 );
    coords.push_back( my_rank + 1 );

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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(),
               hex_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_HEXAHEDRON,
        num_vertices, hex_handle_array, connectivity_array, permutation_list );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
DataTransferKit::MeshContainer<unsigned long int> buildPyramidContainer()
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

    Teuchos::ArrayRCP<unsigned long int> vertex_handle_array(
        vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(),
               vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> pyramid_handle_array(
        pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(),
               pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<unsigned long int> connectivity_array(
        pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(),
               connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return MeshContainer<unsigned long int>(
        vertex_dim, vertex_handle_array, coords_array, DTK_PYRAMID,
        num_vertices, pyramid_handle_array, connectivity_array,
        permutation_list );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshContainer, line_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container = buildLineContainer();

    // Mesh parameters.
    int num_vertices = 2;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }

    // Bounding Boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tri_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container = buildTriContainer();

    // Mesh parameters.
    int num_vertices = 3;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_view[2] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[2] == 1.0 );

    // y
    TEST_ASSERT( coords_view[3] == 0.0 );
    TEST_ASSERT( coords_view[4] == 0.0 );
    TEST_ASSERT( coords_view[5] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[3] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[4] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[5] == 1.0 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }
    // Bounding Boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshContainer, quad_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container = buildQuadContainer();

    // Mesh parameters.
    int num_vertices = 4;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_view[2] == 1.0 );
    TEST_ASSERT( coords_view[3] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[2] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[3] == 0.0 );

    // y
    TEST_ASSERT( coords_view[4] == 0.0 );
    TEST_ASSERT( coords_view[5] == 0.0 );
    TEST_ASSERT( coords_view[6] == 1.0 );
    TEST_ASSERT( coords_view[7] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[4] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[5] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[6] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[7] == 1.0 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }

    // Bounding boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Tet mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tet_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container = buildTetContainer();

    // Mesh parameters.
    int num_vertices = 4;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_view[2] == 1.0 );
    TEST_ASSERT( coords_view[3] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[2] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[3] == 0.0 );

    // y
    TEST_ASSERT( coords_view[4] == 0.0 );
    TEST_ASSERT( coords_view[5] == 0.0 );
    TEST_ASSERT( coords_view[6] == 1.0 );
    TEST_ASSERT( coords_view[7] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[4] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[5] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[6] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[7] == 1.0 );

    // z
    TEST_ASSERT( coords_view[8] == 0.0 );
    TEST_ASSERT( coords_view[9] == 0.0 );
    TEST_ASSERT( coords_view[10] == 0.0 );
    TEST_ASSERT( coords_view[11] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[8] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[9] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[10] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[11] == 1.0 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }

    // Bounding Boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, hex_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container = buildHexContainer();

    // Mesh parameters.
    int num_vertices = 8;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_view[2] == 1.0 );
    TEST_ASSERT( coords_view[3] == 0.0 );
    TEST_ASSERT( coords_view[4] == 0.0 );
    TEST_ASSERT( coords_view[5] == 1.0 );
    TEST_ASSERT( coords_view[6] == 1.0 );
    TEST_ASSERT( coords_view[7] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[2] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[3] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[4] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[5] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[6] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[7] == 0.0 );

    // y
    TEST_ASSERT( coords_view[8] == 0.0 );
    TEST_ASSERT( coords_view[9] == 0.0 );
    TEST_ASSERT( coords_view[10] == 1.0 );
    TEST_ASSERT( coords_view[11] == 1.0 );
    TEST_ASSERT( coords_view[12] == 0.0 );
    TEST_ASSERT( coords_view[13] == 0.0 );
    TEST_ASSERT( coords_view[14] == 1.0 );
    TEST_ASSERT( coords_view[15] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[8] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[9] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[10] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[11] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[12] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[13] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[14] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[15] == 1.0 );

    // z
    TEST_ASSERT( coords_view[16] == 0.0 );
    TEST_ASSERT( coords_view[17] == 0.0 );
    TEST_ASSERT( coords_view[18] == 0.0 );
    TEST_ASSERT( coords_view[19] == 0.0 );
    TEST_ASSERT( coords_view[20] == 1.0 );
    TEST_ASSERT( coords_view[21] == 1.0 );
    TEST_ASSERT( coords_view[22] == 1.0 );
    TEST_ASSERT( coords_view[23] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[16] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[17] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[18] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[19] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[20] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[21] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[22] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[23] == 1.0 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }

    // Bounding Boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
TEUCHOS_UNIT_TEST( MeshContainer, pyramid_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container = buildPyramidContainer();

    // Mesh parameters.
    int num_vertices = 5;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_view[2] == 1.0 );
    TEST_ASSERT( coords_view[3] == 0.0 );
    TEST_ASSERT( coords_view[4] == 0.5 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[2] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[3] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[4] == 0.5 );

    // y
    TEST_ASSERT( coords_view[5] == 0.0 );
    TEST_ASSERT( coords_view[6] == 0.0 );
    TEST_ASSERT( coords_view[7] == 1.0 );
    TEST_ASSERT( coords_view[8] == 1.0 );
    TEST_ASSERT( coords_view[9] == 0.5 );
    TEST_ASSERT( coords_nonconst_view[5] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[6] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[7] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[8] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[9] == 0.5 );

    // z
    TEST_ASSERT( coords_view[10] == 0.0 );
    TEST_ASSERT( coords_view[11] == 0.0 );
    TEST_ASSERT( coords_view[12] == 0.0 );
    TEST_ASSERT( coords_view[13] == 0.0 );
    TEST_ASSERT( coords_view[14] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[10] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[11] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[12] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[13] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[14] == 1.0 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }

    // Bounding Boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
TEUCHOS_UNIT_TEST( MeshContainer, wedge_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container = buildWedgeContainer();

    // Mesh parameters.
    int num_vertices = 6;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_view[2] == 0.5 );
    TEST_ASSERT( coords_view[3] == 0.0 );
    TEST_ASSERT( coords_view[4] == 1.0 );
    TEST_ASSERT( coords_view[5] == 0.5 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[2] == 0.5 );
    TEST_ASSERT( coords_nonconst_view[3] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[4] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[5] == 0.5 );

    // y
    TEST_ASSERT( coords_view[6] == 0.0 );
    TEST_ASSERT( coords_view[7] == 0.0 );
    TEST_ASSERT( coords_view[8] == 1.0 );
    TEST_ASSERT( coords_view[9] == 0.0 );
    TEST_ASSERT( coords_view[10] == 0.0 );
    TEST_ASSERT( coords_view[11] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[6] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[7] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[8] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[9] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[10] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[11] == 1.0 );

    // z
    TEST_ASSERT( coords_view[12] == 0.0 );
    TEST_ASSERT( coords_view[13] == 0.0 );
    TEST_ASSERT( coords_view[14] == 0.0 );
    TEST_ASSERT( coords_view[15] == 1.0 );
    TEST_ASSERT( coords_view[16] == 1.0 );
    TEST_ASSERT( coords_view[17] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[12] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[13] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[14] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[15] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[16] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[17] == 1.0 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }

    // Bounding Boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, parallel_hex_tools_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a mesh container.
    MeshContainer<unsigned long int> mesh_container =
        buildParallelHexContainer();

    // Mesh parameters.
    int num_vertices = 8;

    // Basic container info.
    typedef MeshTools<MeshContainer<unsigned long int>> Tools;
    TEST_ASSERT( Tools::numElements( mesh_container ) == 1 );
    TEST_ASSERT( (int)Tools::numVertices( mesh_container ) == num_vertices );

    // Vertices.
    Teuchos::ArrayRCP<const unsigned long int> vertices_view =
        Tools::verticesView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> vertices_nonconst_view =
        Tools::verticesNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)vertices_view[i] == i );
        TEST_ASSERT( (int)vertices_nonconst_view[i] == i );
    }

    // Coords.
    Teuchos::ArrayRCP<const double> coords_view =
        Tools::coordsView( mesh_container );
    Teuchos::ArrayRCP<double> coords_nonconst_view =
        Tools::coordsNonConstView( mesh_container );
    // x
    TEST_ASSERT( coords_view[0] == 0.0 );
    TEST_ASSERT( coords_view[1] == 1.0 );
    TEST_ASSERT( coords_view[2] == 1.0 );
    TEST_ASSERT( coords_view[3] == 0.0 );
    TEST_ASSERT( coords_view[4] == 0.0 );
    TEST_ASSERT( coords_view[5] == 1.0 );
    TEST_ASSERT( coords_view[6] == 1.0 );
    TEST_ASSERT( coords_view[7] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[0] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[1] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[2] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[3] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[4] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[5] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[6] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[7] == 0.0 );

    // y
    TEST_ASSERT( coords_view[8] == 0.0 );
    TEST_ASSERT( coords_view[9] == 0.0 );
    TEST_ASSERT( coords_view[10] == 1.0 );
    TEST_ASSERT( coords_view[11] == 1.0 );
    TEST_ASSERT( coords_view[12] == 0.0 );
    TEST_ASSERT( coords_view[13] == 0.0 );
    TEST_ASSERT( coords_view[14] == 1.0 );
    TEST_ASSERT( coords_view[15] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[8] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[9] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[10] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[11] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[12] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[13] == 0.0 );
    TEST_ASSERT( coords_nonconst_view[14] == 1.0 );
    TEST_ASSERT( coords_nonconst_view[15] == 1.0 );

    // z
    TEST_ASSERT( coords_view[16] == my_rank );
    TEST_ASSERT( coords_view[17] == my_rank );
    TEST_ASSERT( coords_view[18] == my_rank );
    TEST_ASSERT( coords_view[19] == my_rank );
    TEST_ASSERT( coords_view[20] == my_rank + 1 );
    TEST_ASSERT( coords_view[21] == my_rank + 1 );
    TEST_ASSERT( coords_view[22] == my_rank + 1 );
    TEST_ASSERT( coords_view[23] == my_rank + 1 );
    TEST_ASSERT( coords_nonconst_view[16] == my_rank );
    TEST_ASSERT( coords_nonconst_view[17] == my_rank );
    TEST_ASSERT( coords_nonconst_view[18] == my_rank );
    TEST_ASSERT( coords_nonconst_view[19] == my_rank );
    TEST_ASSERT( coords_nonconst_view[20] == my_rank + 1 );
    TEST_ASSERT( coords_nonconst_view[21] == my_rank + 1 );
    TEST_ASSERT( coords_nonconst_view[22] == my_rank + 1 );
    TEST_ASSERT( coords_nonconst_view[23] == my_rank + 1 );

    // Elements.
    Teuchos::ArrayRCP<const unsigned long int> elements_view =
        Tools::elementsView( mesh_container );
    TEST_ASSERT( elements_view[0] == 12 );
    Teuchos::ArrayRCP<unsigned long int> elements_nonconst_view =
        Tools::elementsNonConstView( mesh_container );
    TEST_ASSERT( elements_nonconst_view[0] == 12 );

    // Connectivity.
    Teuchos::ArrayRCP<const unsigned long int> connectivity_view =
        Tools::connectivityView( mesh_container );
    Teuchos::ArrayRCP<unsigned long int> connectivity_nonconst_view =
        Tools::connectivityNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)connectivity_view[i] == i );
        TEST_ASSERT( (int)connectivity_nonconst_view[i] == i );
    }

    // Permutation.
    Teuchos::ArrayRCP<const int> permutation_view =
        Tools::permutationView( mesh_container );
    Teuchos::ArrayRCP<int> permutation_nonconst_view =
        Tools::permutationNonConstView( mesh_container );
    for ( int i = 0; i < num_vertices; ++i )
    {
        TEST_ASSERT( (int)permutation_view[i] == i );
        TEST_ASSERT( (int)permutation_nonconst_view[i] == i );
    }

    // Bounding Boxes.
    BoundingBox local_box = Tools::localBoundingBox( mesh_container );
    Teuchos::Tuple<double, 6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == my_rank );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == my_rank + 1 );

    BoundingBox global_box =
        Tools::globalBoundingBox( mesh_container, getDefaultComm<int>() );
    Teuchos::Tuple<double, 6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == my_size );
}

//---------------------------------------------------------------------------//
// end tstMeshTools.cpp
//---------------------------------------------------------------------------//
