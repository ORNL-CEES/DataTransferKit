//---------------------------------------------------------------------------//
/*!
 * \file tstTopologyTools.cpp
 * \author Stuart R. Slattery
 * \brief Topology tools unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_TopologyTools.hpp>
#include <DTK_MeshTypes.hpp>

#include <MBInterface.hpp>
#include <MBCore.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
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
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Line point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, line_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    error = moab->set_dimension( 1 );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    double vertex_0[1] = { -1.43 };
    double vertex_1[1] = { 2.98 };

    std::vector<moab::EntityHandle> vertices(2);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle line_segment;
    moab->create_element( moab::MBEDGE, &vertices[0], 2, line_segment );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
    		     moab->type_from_handle( line_segment ) ) == 2 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(1);
    point_0[0] = 0.5;
    Teuchos::Array<double> point_1(1); 
    point_1[0] = 5.2;
    Teuchos::Array<double> point_2(1);
    point_2[0] = 2.98001;
    Teuchos::Array<double> point_3(1);
    point_3[0] = -1.5;
    Teuchos::Array<double> point_4(1);
    point_4[0] = 0.25;
    Teuchos::Array<double> point_5(1);
    point_5[0] = 2.98;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, line_segment, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, line_segment, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, line_segment, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, line_segment, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, line_segment, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, line_segment, moab ) );
}

//---------------------------------------------------------------------------//
// Quad point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, quad_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    error = moab->set_dimension( 2 );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    double vertex_0[2] = { -1.43, -2.3 };
    double vertex_1[2] = { 2.98, -2.12 };
    double vertex_2[2] = { 0.43, 4.2 };
    double vertex_3[2] = { -3.3, 1.1 };

    std::vector<moab::EntityHandle> vertices(3);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_3, vertices[3] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle quadrilateral;
    moab->create_element( moab::MBQUAD, &vertices[0], 4, quadrilateral );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( quadrilateral ) ) == 4 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(2);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    Teuchos::Array<double> point_1(2); 
    point_1[0] = 5.2;
    point_1[1] = 0.9;
    Teuchos::Array<double> point_2(2);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    Teuchos::Array<double> point_3(2);
    point_3[0] = -11.5;
    point_3[1] = 1.5;
    Teuchos::Array<double> point_4(2);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    Teuchos::Array<double> point_5(2);
    point_5[0] = 0.43;
    point_5[1] = 4.2;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, quadrilateral, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, quadrilateral, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, quadrilateral, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, quadrilateral, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, quadrilateral, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, quadrilateral, moab ) );
}

//---------------------------------------------------------------------------//
// Tri point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, tri_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    error = moab->set_dimension( 2 );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    double vertex_0[2] = { -1.43, -2.3 };
    double vertex_1[2] = { 2.98, -2.12 };
    double vertex_2[2] = { 0.43, 4.2 };

    std::vector<moab::EntityHandle> vertices(3);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle triangle;
    moab->create_element( moab::MBTRI, &vertices[0], 3, triangle );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( triangle ) ) == 3 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(2);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    Teuchos::Array<double> point_1(2); 
    point_1[0] = 5.2;
    point_1[1] = 0.9;
    Teuchos::Array<double> point_2(2);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    Teuchos::Array<double> point_3(2);
    point_3[0] = -1.5;
    point_3[1] = 1.5;
    Teuchos::Array<double> point_4(2);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    Teuchos::Array<double> point_5(2);
    point_5[0] = 0.43;
    point_5[1] = 4.2;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, triangle, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, triangle, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, triangle, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, triangle, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, triangle, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, triangle, moab ) );
}

//---------------------------------------------------------------------------//
// Tet point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, tet_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );

    double vertex_0[3] = { -1.43, -2.3, 3.4 };
    double vertex_1[3] = { 2.98, -2.12, 4.3 };
    double vertex_2[3] = { 0.43, 4.2, 2.4 };
    double vertex_3[3] = { 0.98, 3.77, 9.8 };

    std::vector<moab::EntityHandle> vertices(4);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_3, vertices[3] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle tetrahedron;
    moab->create_element( moab::MBTET, &vertices[0], 4, tetrahedron );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( tetrahedron ) ) == 4 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    point_0[2] = 6.2;
    Teuchos::Array<double> point_1(3); 
    point_1[0] = 0.2;
    point_1[1] = 0.9;
    point_1[2] = 1.32;
    Teuchos::Array<double> point_2(3);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    point_2[2] = 9.5;
    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.1;
    point_3[1] = 1.5;
    point_3[2] = -4.8;
    Teuchos::Array<double> point_4(3);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    point_4[2] = 5.4;
    Teuchos::Array<double> point_5(3);
    point_5[0] = 0.98;
    point_5[1] = 3.77;
    point_5[2] = 9.8;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, tetrahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, tetrahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, tetrahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, tetrahedron, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, tetrahedron, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, tetrahedron, moab ) );
}

//---------------------------------------------------------------------------//
// Pyramid point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, pyramid_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );

    double vertex_0[3] = { -1.43, -2.3, 3.4 };
    double vertex_1[3] = { 2.98, -2.12, 4.3 };
    double vertex_2[3] = { 0.43, 4.2, 2.4 };
    double vertex_3[3] = { -3.3, 1.1, 1.1 };
    double vertex_4[3] = { 0.98, 3.77, 9.8 };

    std::vector<moab::EntityHandle> vertices(5);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_3, vertices[3] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_4, vertices[4] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle pyramid;
    moab->create_element( moab::MBPYRAMID, &vertices[0], 5, pyramid );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( pyramid ) ) == 5 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    point_0[2] = 6.2;
    Teuchos::Array<double> point_1(3); 
    point_1[0] = 0.2;
    point_1[1] = 0.9;
    point_1[2] = 1.32;
    Teuchos::Array<double> point_2(3);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    point_2[2] = 9.5;
    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.1;
    point_3[1] = 1.5;
    point_3[2] = -4.8;
    Teuchos::Array<double> point_4(3);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    point_4[2] = 5.2;
    Teuchos::Array<double> point_5(3);
    point_5[0] = 0.98;
    point_5[1] = 3.77;
    point_5[2] = 9.8;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, pyramid, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, pyramid, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, pyramid, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, pyramid, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, pyramid, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, pyramid, moab ) );
}

//---------------------------------------------------------------------------//
// Hexahedron point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, hexahedron_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );

    double vertex_0[3] = { -1.43, -2.3, 3.4 };
    double vertex_1[3] = { 2.98, -2.12, 4.3 };
    double vertex_2[3] = { 0.43, 4.2, 2.4 };
    double vertex_3[3] = { -3.3, 1.1, 1.1 };
    double vertex_4[3] = { -3.43, -2.3, 9.4 };
    double vertex_5[3] = { 1.98, -2.12, 8.3 };
    double vertex_6[3] = { 0.43, 4.2, 7.4 };
    double vertex_7[3] = { -2.3, 1.1, 6.1 };

    std::vector<moab::EntityHandle> vertices(8);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_3, vertices[3] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_4, vertices[4] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_5, vertices[5] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_6, vertices[6] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_7, vertices[7] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle hexahedron;
    moab->create_element( moab::MBHEX, &vertices[0], 8, hexahedron );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( hexahedron ) ) == 8 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    point_0[2] = 5.98;
    Teuchos::Array<double> point_1(3); 
    point_1[0] = 0.2;
    point_1[1] = 0.9;
    point_1[2] = 1.32;
    Teuchos::Array<double> point_2(3);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    point_2[2] = 9.5;
    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.1;
    point_3[1] = 1.5;
    point_3[2] = -4.8;
    Teuchos::Array<double> point_4(3);
    point_4[0] = 1.0;
    point_4[1] = 1.0;
    point_4[2] = 6.0;
    Teuchos::Array<double> point_5(3);
    point_5[0] = -2.3;
    point_5[1] = 1.1;
    point_5[2] = 6.1;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, hexahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, hexahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, hexahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, hexahedron, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, hexahedron, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, hexahedron, moab ) );
}

//---------------------------------------------------------------------------//
// end tstTopologyTools.cpp
//---------------------------------------------------------------------------//

