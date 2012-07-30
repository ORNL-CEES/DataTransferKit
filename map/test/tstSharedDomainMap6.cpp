//---------------------------------------------------------------------------//
/*!
 * \file tstSharedDomainMap6.cpp
 * \author Stuart R. Slattery
 * \brief Shared domain map unit test 6.
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

#include <DTK_SharedDomainMap.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_RendezvousMesh.hpp>

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
    { return DTK_PYRAMID; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 5; }


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
    int num_nodes = edge_length*edge_length*2 + (edge_length-1)*(edge_length-1);
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
	    idx = i + j*edge_length + edge_length*edge_length;
	    node_handles[ idx ] = (long int) num_nodes*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 1.0;
	}
    }
    for ( int j = 0; j < edge_length-1; ++j )
    {
	for ( int i = 0; i < edge_length-1; ++i )
	{
	    idx = i + j*(edge_length-1) + edge_length*edge_length*2;
	    node_handles[ idx ] = (long int) num_nodes*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1) + 0.5;
	    coords[ num_nodes + idx ] = j + 0.5;
	    coords[ 2*num_nodes + idx ] = 0.5;
	}
    }
    
    // Make the pyramids. 
    int num_elements = (edge_length-1)*(edge_length-1)*6;
    Teuchos::Array<long int> pyr_handles( num_elements );
    Teuchos::Array<long int> pyr_connectivity( 5*num_elements );
    int elem_idx, node_idx;
    int v0, v1, v2, v3, v4, v5, v6, v7, v8;
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
	    v4 = node_idx +                   edge_length*edge_length;
	    v5 = node_idx + 1 +               edge_length*edge_length;
	    v6 = node_idx + 1 + edge_length + edge_length*edge_length;
	    v7 = node_idx +     edge_length + edge_length*edge_length;
	    v8 = i + j*(edge_length-1) + edge_length*edge_length*2;

	    // Pyramid 1.
	    elem_idx = i + j*(edge_length-1);
	    pyr_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    pyr_connectivity[elem_idx]                = node_handles[v0];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v1];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v2];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v3];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/6;
	    pyr_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    pyr_connectivity[elem_idx] 	              = node_handles[v1];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v5];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v2];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 3.
	    elem_idx = i + j*(edge_length-1) + 2*num_elements/6;
	    pyr_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    pyr_connectivity[elem_idx] 	              = node_handles[v2];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v6];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v7];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v3];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 4.
	    elem_idx = i + j*(edge_length-1) + 3*num_elements/6;
	    pyr_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    pyr_connectivity[elem_idx]   	      = node_handles[v4];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v0];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v3];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v7];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 5.
	    elem_idx = i + j*(edge_length-1) + 4*num_elements/6;
	    pyr_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    pyr_connectivity[elem_idx]   	      = node_handles[v4];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v5];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v1];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v0];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 6.
	    elem_idx = i + j*(edge_length-1) + 5*num_elements/6;
	    pyr_handles[elem_idx] = num_elements*my_rank + elem_idx;
	    pyr_connectivity[elem_idx]   	      = node_handles[v4];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v7];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v5];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];
	}
    }

    Teuchos::Array<std::size_t> permutation_list( 5 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return MyMesh( node_handles, coords, pyr_handles, pyr_connectivity,
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
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_test6 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<MyMesh> mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP< MeshManager<MyMesh> > source_mesh_manager = Teuchos::rcp( 
	new MeshManager<MyMesh>( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field manager.
    int num_points = 1000;
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( 
			  buildCoordinateField( my_rank, my_size, 
						num_points, edge_size), comm ) );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<MyMesh,MyField> > source_evaluator = 
    	Teuchos::rcp( new MyEvaluator( mesh_blocks[0], comm ) );

    // Create data target. This target is a scalar.
    int field_size = target_coord_manager->field().size() 
		     / target_coord_manager->field().dim();
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( MyField( field_size, 1 ), comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MyMesh,MyField> shared_domain_map( comm );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data.
    int source_rank;
    long int target_dim_size = target_space_manager->field().size() / 
			       target_space_manager->field().dim();
    for ( long int n = 0; n < target_dim_size; ++n )
    {
	source_rank = std::floor(target_coord_manager->field().getData()[n] 
				 / (edge_size-1));
	TEST_ASSERT( source_rank+1 == 
		     target_space_manager->field().getData()[n] );
    }
}

//---------------------------------------------------------------------------//
// end tstSharedDomainMap6.cpp
//---------------------------------------------------------------------------//

