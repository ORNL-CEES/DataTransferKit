//---------------------------------------------------------------------------//
/*!
 * \file tstRendezvous.cpp
 * \author Stuart R. Slattery
 * \brief Rendezvous unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_Rendezvous.hpp>
#include <DTK_RendezvousMesh.hpp>
#include <DTK_BoundingBox.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
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
// Mesh Implementation
//---------------------------------------------------------------------------//

class MyMesh
{
  public:

    typedef int    handle_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<int>& node_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<int>& quad_handles,
	    const Teuchos::Array<int>& quad_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_quad_handles( quad_handles )
	, d_quad_connectivity( quad_connectivity )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    Teuchos::Array<int>::const_iterator nodesBegin() const
    { return d_node_handles.begin(); }

    Teuchos::Array<int>::const_iterator nodesEnd() const
    { return d_node_handles.end(); }

    Teuchos::Array<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    Teuchos::Array<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    Teuchos::Array<int>::const_iterator quadsBegin() const
    { return d_quad_handles.begin(); }

    Teuchos::Array<int>::const_iterator quadsEnd() const
    { return d_quad_handles.end(); }

    Teuchos::Array<int>::const_iterator connectivityBegin() const
    { return d_quad_connectivity.begin(); }

    Teuchos::Array<int>::const_iterator connectivityEnd() const
    { return d_quad_connectivity.end(); }
    

  private:

    Teuchos::Array<int> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<int> d_quad_handles;
    Teuchos::Array<int> d_quad_connectivity;
};

//---------------------------------------------------------------------------//
// DTK Traits Specializations
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

//---------------------------------------------------------------------------//
// Mesh traits specialization for MyMesh
template<>
class MeshTraits<MyMesh>
{
  public:

    typedef MyMesh::handle_type handle_type;
    typedef Teuchos::Array<int>::const_iterator const_node_iterator;
    typedef Teuchos::Array<double>::const_iterator const_coordinate_iterator;
    typedef Teuchos::Array<int>::const_iterator const_element_iterator;
    typedef Teuchos::Array<int>::const_iterator const_connectivity_iterator;

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementType( const MyMesh& mesh )
    { return DTK_FACE; }

    static inline std::size_t elementTopology( const MyMesh& mesh )
    { return DTK_QUADRILATERAL; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 4; }


    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.quadsBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.quadsEnd(); }

    static inline const_connectivity_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Mesh create funciton.
//---------------------------------------------------------------------------//
MyMesh buildMyMesh()
{
    int my_rank = getDefaultComm<int>()->getRank();

    // Make some nodes.
    int num_nodes = 10;
    Teuchos::Array<int> node_handles( num_nodes );
    Teuchos::Array<double> coords( 3*num_nodes );

    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles[i] = (num_nodes / 2)*my_rank + i;
    }
    for ( int i = 0; i < num_nodes / 2; ++i )
    {
	coords[ i ] = my_rank;
	coords[ num_nodes + i ] = i;
	coords[ 2*num_nodes + i ] = 0.0;
    }
    for ( int i = num_nodes / 2; i < num_nodes; ++i )
    {
	coords[ i ] = my_rank + 1;
	coords[ num_nodes + i ] = i - num_nodes/2;
	coords[ 2*num_nodes + i ] = 0.0;
    }
    
    // Make the quads.
    int num_quads = 4;
    Teuchos::Array<int> quad_handles( num_quads );
    Teuchos::Array<int> quad_connectivity( 4*num_quads );
    
    for ( int i = 0; i < num_quads; ++i )
    {
	quad_handles[ i ] = num_quads*my_rank + i;
	quad_connectivity[ i ] = node_handles[i];
	quad_connectivity[ num_quads + i ] = node_handles[num_nodes/2 + i];
	quad_connectivity[ 2*num_quads + i ] = node_handles[num_nodes/2 + i + 1];
	quad_connectivity[ 3*num_quads + i ] = node_handles[i+1];
    }

    return MyMesh( node_handles, coords, quad_handles, quad_connectivity );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

// This test will repartition via RCB the following quad mesh on 4 processes.
/*
          *-------*-------*-------*-------*
          |       |       |       |       |
          |   0   |   1   |   2   |   3   |
          |       |       |       |       |
          *-------*-------*-------*-------*
          |       |       |       |       |
          |   0   |   1   |   2   |   3   |
          |       |       |       |       |
          *-------*-------*-------*-------*
          |       |       |       |       |
          |   0   |   1   |   2   |   3   |
          |       |       |       |       |
          *-------*-------*-------*-------*
          |       |       |       |       |
          |   0   |   1   |   2   |   3   |
          |       |       |       |       |
          *-------*-------*-------*-------*

the result considering node overlap should be:

          *-------*-------*-------*-------*
          |       |       |       |       |
          |   0   |   0   |   2   |   2   |
          |       |       |       |       |
          *-------*-------*-------*-------*
          |       |       |       |       |
          |   0   |   0   |   2   |   2   |
          |       |       |       |       |
          *-------*-------*-------*-------*
          |       |       |       |       |
          |   1   |   1   |   3   |   3   |
          |       |       |       |       |
          *-------*-------*-------*-------*
          |       |       |       |       |
          |   1   |   1   |   3   |   3   |
          |       |       |       |       |
          *-------*-------*-------*-------*

*/
TEUCHOS_UNIT_TEST( Rendezvous, rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    BoundingBox box( -100, -100, -100, 100, 100, 100 );

    // Create a mesh.
    MyMesh my_mesh = buildMyMesh();
    std::string input_file_name;
    for ( int i = 0; i < my_size; ++i )
    {
	if ( my_rank == i )
	{
	    std::stringstream convert;
	    convert << i;
	    input_file_name = "rank_" + convert.str()  + "_in.vtk";
	    createRendezvousMesh( my_mesh )->getMoab()->write_mesh(
		input_file_name.c_str() );
	}
	getDefaultComm<int>()->barrier();
    }

    // Create a rendezvous.
    Rendezvous<MyMesh> rendezvous( getDefaultComm<int>(), box );
    rendezvous.build( my_mesh );

    // Check the mesh.
    getDefaultComm<int>()->barrier();
    std::string output_file_name;
    for ( int i = 0; i < my_size; ++i )
    {
	if ( my_rank == i )
	{
	    std::stringstream convert;
	    convert << i;
	    output_file_name = "rank_" +  convert.str() + "_out.vtk";
	    rendezvous.getMesh()->getMoab()->write_mesh( 
		output_file_name.c_str() );
	}
	getDefaultComm<int>()->barrier();
    }
	
    // Check the partitioning with coordinates.

    // Check the kD-tree with coordinates.
}

//---------------------------------------------------------------------------//
// end tstRendezvous.cpp
//---------------------------------------------------------------------------//
