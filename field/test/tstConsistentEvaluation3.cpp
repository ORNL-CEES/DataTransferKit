//---------------------------------------------------------------------------//
/*!
 * \file tstConsistentEvaluation3.cpp
 * \author Stuart R. Slattery
 * \brief Consistent evaluation unit test 3.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include <DTK_ConsistentEvaluation.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_RendezvousMesh.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
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

    typedef long int    global_ordinal_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<global_ordinal_type>& node_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<global_ordinal_type>& element_handles,
	    const Teuchos::Array<global_ordinal_type>& element_connectivity,
	    const Teuchos::Array<std::size_t>& permutation_list )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_element_handles( element_handles )
	, d_element_connectivity( element_connectivity )
	, d_permutation_list( permutation_list )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    Teuchos::Array<global_ordinal_type>::const_iterator nodesBegin() const
    { return d_node_handles.begin(); }

    Teuchos::Array<global_ordinal_type>::const_iterator nodesEnd() const
    { return d_node_handles.end(); }

    Teuchos::Array<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    Teuchos::Array<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    Teuchos::Array<global_ordinal_type>::const_iterator elementsBegin() const
    { return d_element_handles.begin(); }

    Teuchos::Array<global_ordinal_type>::const_iterator elementsEnd() const
    { return d_element_handles.end(); }

    Teuchos::Array<global_ordinal_type>::const_iterator 
    connectivityBegin() const
    { return d_element_connectivity.begin(); }

    Teuchos::Array<global_ordinal_type>::const_iterator 
    connectivityEnd() const
    { return d_element_connectivity.end(); }
    
    Teuchos::Array<std::size_t>::const_iterator permutationBegin() const
    { return d_permutation_list.begin(); }

    Teuchos::Array<std::size_t>::const_iterator permutationEnd() const
    { return d_permutation_list.end(); }


  private:

    Teuchos::Array<global_ordinal_type> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<global_ordinal_type> d_element_handles;
    Teuchos::Array<global_ordinal_type> d_element_connectivity;
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
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_node_iterator;
    typedef Teuchos::Array<double>::const_iterator 
    const_coordinate_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_element_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_connectivity_iterator;
    typedef Teuchos::Array<std::size_t>::const_iterator 
    const_permutation_iterator;


    static inline std::size_t nodeDim( const MyMesh& mesh )
    { return 3; }

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementTopology( const MyMesh& mesh )
    { return DTK_TETRAHEDRON; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 4; }


    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.elementsBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.elementsEnd(); }

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
    typedef double                     value_type;
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
	    if ( std::find( d_mesh.elementsBegin(),
			    d_mesh.elementsEnd(),
			    elements[n] ) != d_mesh.elementsEnd() )
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
MyMesh buildMyMesh( int my_rank, int my_size, int edge_length )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length*2;
    int node_dim = 3;
    Teuchos::Array<long int> node_handles( num_nodes );
    Teuchos::Array<double> coords( node_dim*num_nodes );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    node_handles[ idx ] = (long int) num_nodes*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + num_nodes / 2;
	    node_handles[ idx ] = (long int) num_nodes*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 1.0;
	}
    }
    
    // Make the tetrahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1)*5;
    Teuchos::Array<long int> tet_handles( num_elements );
    Teuchos::Array<long int> tet_connectivity( 4*num_elements );
    int elem_idx, node_idx;
    int v0, v1, v2, v3, v4, v5, v6, v7;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    // Indices.
	    node_idx = i + j*edge_length;
	    v0 = node_idx;
	    v1 = node_idx + 1;
	    v2 = node_idx + 1 + edge_length;
	    v3 = node_idx +     edge_length;
	    v4 = node_idx +                   num_nodes/2;
	    v5 = node_idx + 1 +               num_nodes/2;
	    v6 = node_idx + 1 + edge_length + num_nodes/2;
	    v7 = node_idx +     edge_length + num_nodes/2; 

	    // Tetrahedron 1.
	    elem_idx = i + j*(edge_length-1);
	    tet_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    tet_connectivity[elem_idx]                = node_handles[v0];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v1];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v3];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v4];

	    // Tetrahedron 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/5;
	    tet_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    tet_connectivity[elem_idx] 	              = node_handles[v1];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v2];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v3];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v6];

	    // Tetrahedron 3.
	    elem_idx = i + j*(edge_length-1) + 2*num_elements/5;
	    tet_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    tet_connectivity[elem_idx] 	              = node_handles[v6];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v5];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v4];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v1];

	    // Tetrahedron 4.
	    elem_idx = i + j*(edge_length-1) + 3*num_elements/5;
	    tet_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    tet_connectivity[elem_idx]   	      = node_handles[v4];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v7];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v3];

	    // Tetrahedron 5.
	    elem_idx = i + j*(edge_length-1) + 4*num_elements/5;
	    tet_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    tet_connectivity[elem_idx] 	              = node_handles[v3];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v1];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v4];
	}
    }

    Teuchos::Array<std::size_t> permutation_list( 4 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return MyMesh( node_handles, coords, tet_handles, tet_connectivity,
		   permutation_list );
}

//---------------------------------------------------------------------------//
// Coordinate field create function.
//---------------------------------------------------------------------------//
MyField buildCoordinateField( int my_rank, int my_size, 
			      int num_points, int edge_size )
{
    std::srand( my_rank*num_points*2 );
    int point_dim = 3;
    MyField coordinate_field( num_points*point_dim, point_dim );

    for ( int i = 0; i < num_points; ++i )
    {
	*(coordinate_field.begin() + i) = 
	    my_size * (edge_size-1) * (double) std::rand() / RAND_MAX;
	*(coordinate_field.begin() + num_points + i ) = 
	    (edge_size-1) * (double) std::rand() / RAND_MAX;
	*(coordinate_field.begin() + 2*num_points + i ) = 0.5;
    }

    return coordinate_field;
}

//---------------------------------------------------------------------------//
// Unit test.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentEvaluation, consistent_evaluation_test3 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh.
    int edge_size = 4;
    MyMesh source_mesh = buildMyMesh( my_rank, my_size, edge_size );

    // Setup target coordinate field.
    int num_points = (edge_size-1)*(edge_size-1);
    MyField target_coords = buildCoordinateField( my_rank, my_size, 
						  num_points, edge_size );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<MyMesh,MyField> > my_evaluator = 
    	Teuchos::rcp( new MyEvaluator( source_mesh, comm ) );

    // Create data target.
    MyField::size_type target_size = 
	target_coords.size() / target_coords.dim();
    MyField my_target( target_size, 1 );

    // Setup and apply the consistent evaluation mapping.
    typedef ConsistentEvaluation<MyMesh,MyField> MapType;
    Teuchos::RCP<MapType> consistent_evaluation = 
    	Teuchos::rcp( new MapType( comm ) );
    consistent_evaluation->setup( source_mesh, target_coords );
    consistent_evaluation->apply( my_evaluator, my_target );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data.
    int source_rank;
    for ( long int n = 0; n < my_target.size(); ++n )
    {
	source_rank = std::floor(target_coords.getData()[n] / (edge_size-1));
	TEST_ASSERT( source_rank+1 == my_target.getData()[n] );
    }
}

//---------------------------------------------------------------------------//
// end tstConsistentEvaluation3.cpp
//---------------------------------------------------------------------------//

