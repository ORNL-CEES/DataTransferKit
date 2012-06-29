//---------------------------------------------------------------------------//
/*!
 * \file tstTransferOperator.cpp
 * \author Stuart R. Slattery
 * \brief Transfer operator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_TransferOperator.hpp>
#include <DTK_ConsistentEvaluation.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
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

    typedef int    global_ordinal_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<int>& node_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<int>& quad_handles,
	    const Teuchos::Array<int>& quad_connectivity,
	    const Teuchos::Array<std::size_t>& permutation_list )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_quad_handles( quad_handles )
	, d_quad_connectivity( quad_connectivity )
	, d_permutation_list( permutation_list )
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
    
    Teuchos::Array<std::size_t>::const_iterator permutationBegin() const
    { return d_permutation_list.begin(); }

    Teuchos::Array<std::size_t>::const_iterator permutationEnd() const
    { return d_permutation_list.end(); }


  private:

    Teuchos::Array<int> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<int> d_quad_handles;
    Teuchos::Array<int> d_quad_connectivity;
    Teuchos::Array<std::size_t> d_permutation_list;
};

//---------------------------------------------------------------------------//
// Field implementation.
//---------------------------------------------------------------------------//
class MyField
{
  public:

    typedef double value_type;
    typedef Teuchos::Array<double>::size_type size_type;
    typedef Teuchos::Array<double>::iterator iterator;
    typedef Teuchos::Array<double>::const_iterator const_iterator;

    MyField( size_type size, std::size_t dim )
	: d_dim( dim )
	, d_data( size )
    { /* ... */ }

    ~MyField()
    { /* ... */ }

    std::size_t dim() const
    { return d_dim; }

    size_type size() const
    { return d_data.size(); }

    bool empty() const
    { return d_data.empty(); }

    iterator begin()
    { return d_data.begin(); }

    const_iterator begin() const
    { return d_data.begin(); }

    iterator end()
    { return d_data.end(); }

    const_iterator end() const
    { return d_data.end(); }

    Teuchos::Array<double>& getData()
    { return d_data; }

    const Teuchos::Array<double>& getData() const
    { return d_data; }

  private:
    std::size_t d_dim;
    Teuchos::Array<double> d_data;
};

//---------------------------------------------------------------------------//
// DTK implementations.
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

//---------------------------------------------------------------------------//
// Mesh traits specialization for MyMesh
template<>
class MeshTraits<MyMesh>
{
  public:

    typedef MyMesh::global_ordinal_type global_ordinal_type;
    typedef Teuchos::Array<int>::const_iterator const_node_iterator;
    typedef Teuchos::Array<double>::const_iterator const_coordinate_iterator;
    typedef Teuchos::Array<int>::const_iterator const_element_iterator;
    typedef Teuchos::Array<int>::const_iterator const_connectivity_iterator;
    typedef Teuchos::Array<std::size_t>::const_iterator 
    const_permutation_iterator;

    static inline std::size_t nodeDim( const MyMesh& mesh )
    { return 2; }

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


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

    static inline const_permutation_iterator permutationBegin( const MyMesh& mesh )
    { return mesh.permutationBegin(); }

    static inline const_permutation_iterator permutationEnd( const MyMesh& mesh )
    { return mesh.permutationEnd(); }
};

//---------------------------------------------------------------------------//
// Field Traits specification for MyField
template<>
class FieldTraits<MyField>
{
  public:

    typedef MyField                    field_type;
    typedef double                                    value_type;
    typedef MyField::size_type         size_type;
    typedef MyField::iterator          iterator;
    typedef MyField::const_iterator    const_iterator;

    static inline size_type dim( const MyField& field )
    { return field.dim(); }

    static inline size_type size( const MyField& field )
    { return field.size(); }

    static inline bool empty( const MyField& field )
    { return field.empty(); }

    static inline iterator begin( MyField& field )
    { return field.begin(); }

    static inline const_iterator begin( const MyField& field )
    { return field.begin(); }

    static inline iterator end( MyField& field )
    { return field.end(); }

    static inline const_iterator end( const MyField& field )
    { return field.end(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// FieldEvaluator Implementation.
class MyEvaluator : public DataTransferKit::FieldEvaluator<MyMesh,MyField>
{
  public:

    MyEvaluator( const MyMesh& mesh, 
		 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
	: d_mesh( mesh )
	, d_comm( comm )
    { /* ... */ }

    ~MyEvaluator()
    { /* ... */ }

    MyField evaluate( 
	const Teuchos::ArrayRCP<MyMesh::global_ordinal_type>& elements,
	const Teuchos::ArrayRCP<double>& coords )
    {
	MyField evaluated_data( elements.size(), 1 );
	for ( int n = 0; n < elements.size(); ++n )
	{
	    if ( std::find( d_mesh.quadsBegin(),
			    d_mesh.quadsEnd(),
			    elements[n] ) != d_mesh.quadsEnd() )
	    {
		*(evaluated_data.begin() + n ) = d_comm->getRank() + 1.0;
	    }
	    else
	    {
		*(evaluated_data.begin() + n ) = 0.0;
	    }
	}
	return evaluated_data;
    }

  private:

    MyMesh d_mesh;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// Mesh create function.
//---------------------------------------------------------------------------//
/*
  Make the following mesh partitioned on 4 processors.

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

 */
MyMesh buildMyMesh()
{
    int my_rank = getDefaultComm<int>()->getRank();

    // Make some nodes.
    int num_nodes = 10;
    int node_dim = 2;
    Teuchos::Array<int> node_handles( num_nodes );
    Teuchos::Array<double> coords( node_dim*num_nodes );

    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles[i] = (num_nodes / 2)*my_rank + i;
    }
    for ( int i = 0; i < num_nodes / 2; ++i )
    {
	coords[ i ] = my_rank;
	coords[ num_nodes + i ] = i;
    }
    for ( int i = num_nodes / 2; i < num_nodes; ++i )
    {
	coords[ i ] = my_rank + 1;
	coords[ num_nodes + i ] = i - num_nodes/2;
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

    Teuchos::Array<std::size_t> permutation_list( 4 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return MyMesh( node_handles, coords, quad_handles, quad_connectivity,
		   permutation_list );
}

//---------------------------------------------------------------------------//
// Coordinate field create function.
//---------------------------------------------------------------------------//
/*
  Make the following coordinate field partitioned on 4 processors.

   ------- ------- ------- -------
  |       |       |       |       |
  |   *   |   *   |   *   |   *   |
  |   3   |   3   |   3   |   3   |
   ------- ------- ------- ------- 
  |       |       |       |       |
  |   *   |   *   |   *   |   *   |
  |   2   |   2   |   2   |   2   |
   ------- ------- ------- ------- 
  |       |       |       |       |
  |   *   |   *   |   *   |   *   |
  |   1   |   1   |   1   |   1   |
   ------- ------- ------- ------- 
  |       |       |       |       |
  |   *   |   *   |   *   |   *   |
  |   0   |   0   |   0   |   0   |
   ------- ------- ------- ------- 

 */
MyField buildCoordinateField()
{
    int num_points = 4;
    int point_dim = 2;
    MyField coordinate_field( num_points*point_dim, point_dim );

    for ( int i = 0; i < num_points; ++i )
    {
	*(coordinate_field.begin() + i) = i + 0.5;
	*(coordinate_field.begin() + num_points + i ) = 
	    getDefaultComm<int>()->getRank() + 0.5;
    }

    return coordinate_field;
}

//---------------------------------------------------------------------------//
// Unit tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( TransferOperator, transfer_operator_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_size = comm->getSize();

    // This is a 4 processor test.
    if ( my_size == 4 )
    {

	// Setup source mesh.
	Teuchos::ArrayRCP<MyMesh> mesh_blocks( 1 );
	mesh_blocks[0] = buildMyMesh();
	Teuchos::RCP< MeshManager<MyMesh> > mesh_manager = Teuchos::rcp( 
	    new MeshManager<MyMesh>( mesh_blocks, comm, 2 ) );

	// Setup target coordinate field.
	MyField target_coords = buildCoordinateField();

	// Create field evaluator.
	Teuchos::RCP< FieldEvaluator<MyMesh,MyField> > my_evaluator = 
	    Teuchos::rcp( new MyEvaluator( mesh_blocks[0], comm ) );

	// Create data target.
	int target_size = target_coords.size() / target_coords.dim();
	MyField my_target( target_size, 1 );

	// Setup consistent evaluation mapping.
	typedef ConsistentEvaluation<MyMesh,MyField> MapType;
	Teuchos::RCP<MapType> consistent_evaluation = 
	    Teuchos::rcp( new MapType( comm ) );

	// Setup and apply the transfer operator to the fields.
	TransferOperator<MapType> transfer_operator( consistent_evaluation );
	transfer_operator.setup( mesh_manager, target_coords );
	transfer_operator.apply( my_evaluator, my_target );

	// Check the data transfer.
	for ( int n = 0; n < my_target.size(); ++n )
	{
	    TEST_ASSERT( *(my_target.begin()+n) == n + 1 );
	}

	// Clear the data.
	my_target.getData().clear();
	my_target.getData().resize( target_size );
	for ( int n = 0; n < my_target.size(); ++n )
	{
	    TEST_ASSERT( *(my_target.begin()+n) == 0.0 );
	}

	// Apply the transfer again.
	transfer_operator.apply( my_evaluator, my_target );

	// Check the data transfer again.
	for ( int n = 0; n < my_target.size(); ++n )
	{
	    TEST_ASSERT( *(my_target.begin()+n) == n + 1 );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//

