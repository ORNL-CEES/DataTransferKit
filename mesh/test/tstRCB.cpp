//---------------------------------------------------------------------------//
/*!
 * \file tstRCB.cpp
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

#include <DTK_RCB.hpp>
#include <DTK_BoundingBox.hpp>
#include <DTK_CoreTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_Exception.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
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
// Mesh Implementation
//---------------------------------------------------------------------------//

class MyMesh
{
  public:

    typedef int    handle_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const std::vector<int>& node_handles,
	    const std::vector<double>& coords,
	    const std::vector<int>& hex_handles,
	    const std::vector<int>& hex_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_hex_handles( hex_handles )
	, d_hex_connectivity( hex_connectivity )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    std::vector<int>::const_iterator nodesBegin() const
    { return d_node_handles.begin(); }

    std::vector<int>::const_iterator nodesEnd() const
    { return d_node_handles.end(); }

    std::vector<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    std::vector<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    std::vector<int>::const_iterator hexesBegin() const
    { return d_hex_handles.begin(); }

    std::vector<int>::const_iterator hexesEnd() const
    { return d_hex_handles.end(); }

    std::vector<int>::const_iterator connectivityBegin() const
    { return d_hex_connectivity.begin(); }

    std::vector<int>::const_iterator connectivityEnd() const
    { return d_hex_connectivity.end(); }
    

  private:

    std::vector<int> d_node_handles;
    std::vector<double> d_coords;
    std::vector<int> d_hex_handles;
    std::vector<int> d_hex_connectivity;
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
    typedef std::vector<int>::const_iterator const_node_iterator;
    typedef std::vector<double>::const_iterator const_coordinate_iterator;
    typedef std::vector<int>::const_iterator const_element_iterator;
    typedef std::vector<int>::const_iterator const_connectivity_iterator;
    

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementType( const MyMesh& mesh )
    { return DTK_REGION; }

    static inline std::size_t elementTopology( const MyMesh& mesh )
    { return DTK_HEXAHEDRON; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 8; }

    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.hexesBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.hexesEnd(); }

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
    int my_size = getDefaultComm<int>()->getSize();

    // Make some random nodes.
    std::srand( 1 );
    std::vector<double> random_numbers;
    int num_rand = 3*my_size;
    for ( int i = 0; i < num_rand; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    std::vector<int> node_handles;
    std::vector<double> coords;
    for ( int i = 0; i < my_size; ++i )
    {
	node_handles.push_back( i*my_size + my_rank );
	coords.push_back( random_numbers[3*i] );
	coords.push_back( random_numbers[3*i+1] );
	coords.push_back( random_numbers[3*i+2] );
    }

    // Empty element vectors. We only need nodes for these tests.
    std::vector<int> element_handles;
    std::vector<int> element_connectivity;
    
    return MyMesh( node_handles, coords, element_handles, element_connectivity );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( RCB, rcb_test )
{
    using namespace DataTransferKit;

    // Create my mesh.
    typedef MeshTraits<MyMesh> MT;
    MyMesh my_mesh = buildMyMesh();

    // All of the nodes will be partitioned.
    int num_nodes = 3 * std::distance( MT::nodesBegin( my_mesh ),
				      MT::nodesEnd( my_mesh ) );
    std::vector<char> active_nodes( num_nodes, 1 );

    // Partition the mesh with RCB.
    typedef RCB<MyMesh>::zoltan_id_type zoltan_id_type;
    RCB<MyMesh> rcb( my_mesh, active_nodes, getDefaultComm<int>() );
    rcb.partition();

    // Get the random numbers that were used to compute the node coordinates.
    std::srand( 1 );
    std::vector<double> random_numbers;
    for ( int i = 0; i < num_nodes; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the nodes.
    typename MT::const_coordinate_iterator coord_iterator;
    std::vector<double>::const_iterator rand_iterator = random_numbers.begin();
    for ( coord_iterator = MT::coordsBegin( my_mesh );
	  coord_iterator != MT::coordsEnd( my_mesh );
	  ++coord_iterator, ++rand_iterator )
    {
	TEST_ASSERT( *coord_iterator == *rand_iterator );
    }

    // Get MPI parameters.
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Check import parameters.
    int num_import = rcb.getNumImport();

    Teuchos::ArrayView<zoltan_id_type> import_global_ids = 
	rcb.getImportGlobalIds();
    TEST_ASSERT( import_global_ids.size() == num_import );

    Teuchos::ArrayView<zoltan_id_type> import_local_ids = 
	rcb.getImportLocalIds();
    TEST_ASSERT( import_local_ids.size() == num_import );

    Teuchos::ArrayView<int> import_procs = rcb.getImportProcs();
    TEST_ASSERT( import_procs.size() == num_import );

    Teuchos::ArrayView<int> import_parts = rcb.getImportParts();
    TEST_ASSERT( import_parts.size() == num_import );

    for ( int i = 0; i < num_import; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( import_procs[i] != my_rank &&
		     import_procs[i] >= 0 &&
		     import_procs[i] < my_size );

	TEST_ASSERT( import_parts[i] == my_rank );
    }

    // Check export parameters.
    int num_export = rcb.getNumExport();

    Teuchos::ArrayView<zoltan_id_type> export_global_ids = 
	rcb.getExportGlobalIds();
    TEST_ASSERT( export_global_ids.size() == num_export );

    Teuchos::ArrayView<zoltan_id_type> export_local_ids = 
	rcb.getExportLocalIds();
    TEST_ASSERT( export_local_ids.size() == num_export );

    Teuchos::ArrayView<int> export_procs = rcb.getExportProcs();
    TEST_ASSERT( export_procs.size() == num_export );
    
    Teuchos::ArrayView<int> export_parts = rcb.getExportParts();
    TEST_ASSERT( export_parts.size() == num_export );

    for ( int i = 0; i < num_export; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( export_procs[i] == export_parts[i] );

	TEST_ASSERT( export_procs[i] != my_rank &&
		     export_procs[i] >= 0 &&
		     export_procs[i] < my_size );

	TEST_ASSERT( export_parts[i] != my_rank &&
		     export_parts[i] >= 0 &&
		     export_parts[i] < my_size );
    }
}

//---------------------------------------------------------------------------//
// end tstRCB.cpp
//---------------------------------------------------------------------------//
