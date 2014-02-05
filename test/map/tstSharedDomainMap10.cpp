//---------------------------------------------------------------------------//
/*!
 * \file tstSharedDomainMap10.cpp
 * \author Stuart R. Slattery
 * \brief Shared domain map unit test 10 for 2D hybrid meshes.
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
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
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
// Field implementation.
//---------------------------------------------------------------------------//
class MyField
{
  public:

    typedef double value_type;
    typedef Teuchos::Array<double>::size_type size_type;
    typedef Teuchos::Array<double>::iterator iterator;
    typedef Teuchos::Array<double>::const_iterator const_iterator;

    MyField( size_type size, int dim )
	: d_dim( dim )
	, d_data( dim*size, 0.0 )
    { /* ... */ }

    ~MyField()
    { /* ... */ }

    int dim() const
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
    int d_dim;
    Teuchos::Array<double> d_data;
};

//---------------------------------------------------------------------------//
// DTK implementations.
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

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
class MyEvaluator : 
    public DataTransferKit::FieldEvaluator<int,MyField>
{
  public:

    MyEvaluator( const DataTransferKit::MeshContainer<int>& mesh, 
		 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
	: d_mesh( mesh )
	, d_comm( comm )
    { /* ... */ }

    ~MyEvaluator()
    { /* ... */ }

    MyField evaluate( 
	const Teuchos::ArrayRCP<
	DataTransferKit::MeshContainer<int>::global_ordinal_type>& elements,
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

    DataTransferKit::MeshContainer<int>  d_mesh;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// Mesh create functions.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> >
buildTriMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some vertices.
    int num_vertices = edge_length*edge_length;
    int vertex_dim = 2;
    Teuchos::ArrayRCP<int> vertex_handles( num_vertices );
    Teuchos::ArrayRCP<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    vertex_handles[ idx ] = (int) num_vertices*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j;
	}
    }
    
    // Make the triangles. 
    int num_elements = (edge_length-1)*(edge_length-1)*2;
    Teuchos::ArrayRCP<int> tri_handles( num_elements );
    Teuchos::ArrayRCP<int> tri_connectivity( 3*num_elements );
    int elem_idx, vertex_idx;
    int v0, v1, v2, v3;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    // Indices.
	    vertex_idx = i + j*edge_length;
	    v0 = vertex_idx;
	    v1 = vertex_idx + 1;
	    v2 = vertex_idx + 1 + edge_length;
	    v3 = vertex_idx +     edge_length;

	    // Triangle 1.
	    elem_idx = i + j*(edge_length-1);
	    tri_handles[elem_idx] = elem_idx + elem_offset;
	    tri_connectivity[elem_idx]                = vertex_handles[v0];
	    tri_connectivity[num_elements+elem_idx]   = vertex_handles[v1];
	    tri_connectivity[2*num_elements+elem_idx] = vertex_handles[v2];

	    // Triangle 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/2;
	    tri_handles[elem_idx] = elem_idx + elem_offset;
	    tri_connectivity[elem_idx] 	              = vertex_handles[v2];
	    tri_connectivity[num_elements+elem_idx]   = vertex_handles[v3];
	    tri_connectivity[2*num_elements+elem_idx] = vertex_handles[v0];
	}
    }

    Teuchos::ArrayRCP<int> permutation_list( 3 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp( 
	new DataTransferKit::MeshContainer<int>( 2, vertex_handles, coords, 
						 DataTransferKit::DTK_TRIANGLE, 3,
						 tri_handles, tri_connectivity,
						 permutation_list ) );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> >
buildTiledTriMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some vertices.
    int num_vertices = edge_length*edge_length;
    int vertex_dim = 2;
    Teuchos::ArrayRCP<int> vertex_handles( num_vertices );
    Teuchos::ArrayRCP<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    vertex_handles[ idx ] = (int) num_vertices*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j + my_rank*(edge_length-1);
	}
    }
    
    // Make the triangles. 
    int num_elements = (edge_length-1)*(edge_length-1)*2;
    Teuchos::ArrayRCP<int> tri_handles( num_elements );
    Teuchos::ArrayRCP<int> tri_connectivity( 3*num_elements );
    int elem_idx, vertex_idx;
    int v0, v1, v2, v3;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    // Indices.
	    vertex_idx = i + j*edge_length;
	    v0 = vertex_idx;
	    v1 = vertex_idx + 1;
	    v2 = vertex_idx + 1 + edge_length;
	    v3 = vertex_idx +     edge_length;

	    // Triangle 1.
	    elem_idx = i + j*(edge_length-1);
	    tri_handles[elem_idx] = elem_idx + elem_offset;
	    tri_connectivity[elem_idx]                = vertex_handles[v0];
	    tri_connectivity[num_elements+elem_idx]   = vertex_handles[v1];
	    tri_connectivity[2*num_elements+elem_idx] = vertex_handles[v2];

	    // Triangle 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/2;
	    tri_handles[elem_idx] = elem_idx + elem_offset;
	    tri_connectivity[elem_idx] 	              = vertex_handles[v2];
	    tri_connectivity[num_elements+elem_idx]   = vertex_handles[v3];
	    tri_connectivity[2*num_elements+elem_idx] = vertex_handles[v0];
	}
    }

    Teuchos::ArrayRCP<int> permutation_list( 3 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp(
	new DataTransferKit::MeshContainer<int>( 2, vertex_handles, coords, 
						 DataTransferKit::DTK_TRIANGLE, 3,
						 tri_handles, tri_connectivity,
						 permutation_list ) );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildNullTriMesh()
{
    Teuchos::ArrayRCP<int> vertex_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> tri_handles(0);
    Teuchos::ArrayRCP<int> tri_connectivity(0);
    Teuchos::ArrayRCP<int> permutation_list(3);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp( 
	new DataTransferKit::MeshContainer<int>( 2, vertex_handles, coords, 
						 DataTransferKit::DTK_TRIANGLE, 3,
						 tri_handles, tri_connectivity,
						 permutation_list ) );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> >  
buildQuadMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some vertices.
    int num_vertices = edge_length*edge_length;
    int vertex_dim = 2;
    Teuchos::ArrayRCP<int> vertex_handles( num_vertices );
    Teuchos::ArrayRCP<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    vertex_handles[ idx ] = (int) num_vertices*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j;
	}
    }
    
    // Make the quadrilaterals. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::ArrayRCP<int> quad_handles( num_elements );
    Teuchos::ArrayRCP<int> quad_connectivity( 4*num_elements );
    int elem_idx, vertex_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    vertex_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    quad_handles[elem_idx] = elem_idx + elem_offset;

	    quad_connectivity[elem_idx] 
		= vertex_handles[vertex_idx];

	    quad_connectivity[num_elements+elem_idx] 
		= vertex_handles[vertex_idx+1];

	    quad_connectivity[2*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length+1];

	    quad_connectivity[3*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length];
	}
    }

    Teuchos::ArrayRCP<int> permutation_list( 4 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp( 
	new DataTransferKit::MeshContainer<int>( 2, vertex_handles, coords, 
						 DataTransferKit::DTK_QUADRILATERAL, 4,
						 quad_handles, quad_connectivity,
						 permutation_list ) );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> >  
buildTiledQuadMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some vertices.
    int num_vertices = edge_length*edge_length;
    int vertex_dim = 2;
    Teuchos::ArrayRCP<int> vertex_handles( num_vertices );
    Teuchos::ArrayRCP<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    vertex_handles[ idx ] = (int) num_vertices*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j + my_rank*(edge_length-1);
	}
    }
    
    // Make the quadrilaterals. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::ArrayRCP<int> quad_handles( num_elements );
    Teuchos::ArrayRCP<int> quad_connectivity( 4*num_elements );
    int elem_idx, vertex_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    vertex_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    quad_handles[elem_idx] = elem_idx + elem_offset;

	    quad_connectivity[elem_idx] 
		= vertex_handles[vertex_idx];

	    quad_connectivity[num_elements+elem_idx] 
		= vertex_handles[vertex_idx+1];

	    quad_connectivity[2*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length+1];

	    quad_connectivity[3*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length];
	}
    }

    Teuchos::ArrayRCP<int> permutation_list( 4 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp( 
	new DataTransferKit::MeshContainer<int>( 2, vertex_handles, coords, 
						 DataTransferKit::DTK_QUADRILATERAL, 4,
						 quad_handles, quad_connectivity,
						 permutation_list ) );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> > 
buildNullQuadMesh()
{
    Teuchos::ArrayRCP<int> vertex_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> quad_handles(0);
    Teuchos::ArrayRCP<int> quad_connectivity(0);
    Teuchos::ArrayRCP<int> permutation_list(4);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp( 
	new DataTransferKit::MeshContainer<int>( 2, vertex_handles, coords, 
						 DataTransferKit::DTK_QUADRILATERAL, 4,
						 quad_handles, quad_connectivity,
						 permutation_list ) );
}

//---------------------------------------------------------------------------//
// Coordinate field create functions.
//---------------------------------------------------------------------------//
void buildCoordinateField( int my_rank, int my_size, 
			   int num_points, int edge_size,
			   Teuchos::RCP<MyField>& coordinate_field )
{
    std::srand( my_rank*num_points*2 );
    for ( int i = 0; i < num_points; ++i )
    {
	*(coordinate_field->begin() + i) = 
	    my_size * (edge_size-1) * (double) std::rand() / RAND_MAX;
	*(coordinate_field->begin() + num_points + i ) = 
	    (edge_size-1) * (double) std::rand() / RAND_MAX;
    }
}

//---------------------------------------------------------------------------//
void buildExpandedCoordinateField( int my_rank, int my_size, 
				   int num_points, int edge_size,
				   Teuchos::RCP<MyField>& coordinate_field )
{
    std::srand( my_rank*num_points*2 );
    for ( int i = 0; i < num_points; ++i )
    {
	*(coordinate_field->begin() + i) = 
	    my_size * (edge_size) * (double) std::rand() / RAND_MAX - 0.5;
	*(coordinate_field->begin() + num_points + i ) = 
	    (edge_size) * (double) std::rand() / RAND_MAX - 0.5;
    }
}

//---------------------------------------------------------------------------//
void buildTiledCoordinateField( int my_rank, int my_size, 
				int num_points, int edge_size,
				Teuchos::RCP<MyField>& coordinate_field )
{
    std::srand( my_rank*num_points*2 );
    for ( int i = 0; i < num_points; ++i )
    {
	*(coordinate_field->begin() + i) = 
	    my_size * (edge_size) * (double) std::rand() / RAND_MAX - 0.5;
	*(coordinate_field->begin() + num_points + i ) = 
	    my_size * (edge_size) * (double) std::rand() / RAND_MAX - 0.5;
    }
}


//---------------------------------------------------------------------------//
// Global test parameters.
//---------------------------------------------------------------------------//

// number of random points to be generated.
int num_points = 1000;

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
// All points will be in the mesh in this test.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_test10 )
{
    using namespace DataTransferKit;
    typedef MeshContainer<int> MeshType;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Compute element ordinal offsets so we make unique global ordinals.
    int edge_size = 10;
    int tri_offset_1 = 0;
    int quad_offset_1 = tri_offset_1 + (edge_size+1)*(edge_size+1)*2;
    int tri_offset_2 = quad_offset_1 + (edge_size+1)*(edge_size+1);
    int quad_offset_2 = tri_offset_2 + (edge_size+1)*(edge_size+1)*2;

    // Setup source mesh manager.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    if ( my_rank == 0 )
    {
	mesh_blocks[0] = 
	    buildTriMesh( my_rank, my_size, edge_size, tri_offset_1 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else if ( my_rank == 1 )
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildQuadMesh( my_rank, my_size, edge_size, quad_offset_1 );
    }
    else if ( my_rank == 2 )
    {
	mesh_blocks[0] = 
	    buildTriMesh( my_rank, my_size, edge_size, tri_offset_2 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else 
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildQuadMesh( my_rank, my_size, edge_size, quad_offset_2 );
    }
    comm->barrier();

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > source_mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 2 ) );

    // Setup target coordinate field manager.
    int point_dim = 2;
    Teuchos::RCP<MyField> coordinate_field = 
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildCoordinateField( my_rank, my_size, num_points, edge_size,
			  coordinate_field );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<int,MyField> > source_evaluator;
    if ( my_rank == 0 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );
    }
    else if ( my_rank == 1 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[1], comm ) );
    }
    else if ( my_rank == 2 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );
    }
    else
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[1], comm ) );
    }
    comm->barrier();

    // Create data target. This target is a scalar.
    int target_dim = 1;
    Teuchos::RCP<MyField> target_field =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MeshType,MyField> shared_domain_map( 
	comm, source_mesh_manager->dim() );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data.
    int source_rank;
    for ( int n = 0; n < num_points; ++n )
    {
	source_rank = std::floor(*(coordinate_field->begin()+n) / (edge_size-1));
	for ( int d = 0; d < target_dim; ++d )
	{
	    TEST_ASSERT( source_rank+1 == 
			 *(target_space_manager->field()->begin()
			   +n+d*num_points) );
	}
    }
}

//---------------------------------------------------------------------------//
// Some points will be outside of the mesh in this test.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_expanded_test10 )
{
    using namespace DataTransferKit;
    typedef MeshContainer<int> MeshType;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Compute element ordinal offsets so we make unique global ordinals.
    int edge_size = 10;
    int tri_offset_1 = 0;
    int quad_offset_1 = tri_offset_1 + (edge_size+1)*(edge_size+1)*2;
    int tri_offset_2 = quad_offset_1 + (edge_size+1)*(edge_size+1);
    int quad_offset_2 = tri_offset_2 + (edge_size+1)*(edge_size+1)*2;

    // Setup source mesh manager.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    if ( my_rank == 0 )
    {
	mesh_blocks[0] = 
	    buildTriMesh( my_rank, my_size, edge_size, tri_offset_1 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else if ( my_rank == 1 )
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildQuadMesh( my_rank, my_size, edge_size, quad_offset_1 );
    }
    else if ( my_rank == 2 )
    {
	mesh_blocks[0] = 
	    buildTriMesh( my_rank, my_size, edge_size, tri_offset_2 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else 
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildQuadMesh( my_rank, my_size, edge_size, quad_offset_2 );
    }
    comm->barrier();

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > source_mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 2 ) );

    // Setup target coordinate field manager.
    int point_dim = 2;
    Teuchos::RCP<MyField> coordinate_field = 
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
				  coordinate_field );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<int,MyField> > source_evaluator;
    if ( my_rank == 0 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );
    }
    else if ( my_rank == 1 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[1], comm ) );
    }
    else if ( my_rank == 2 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );
    }
    else
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[1], comm ) );
    }
    comm->barrier();

    // Create data target. This target is a scalar.
    int target_dim = 1;
    Teuchos::RCP<MyField> target_field =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MeshType ,MyField> shared_domain_map( 
	comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data if it is in the mesh and 0.0 if it is outside.
    int source_rank;
    Teuchos::Array<int> missing_points;
    for ( int n = 0; n < num_points; ++n )
    {
	if ( *(coordinate_field->begin()+n) < 0.0 ||
	     *(coordinate_field->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field->begin()+n+num_points) 
	     > edge_size-1 )
	{
	    missing_points.push_back(n);	
	    for ( int d = 0; d < target_dim; ++d )
	    {
		TEST_ASSERT( 0.0 == *(target_space_manager->field()->begin()
				      +n+d*num_points) );
	    }
	}
	else
	{
	    source_rank = std::floor(target_coord_manager->field()->getData()[n] 
	    			     / (edge_size-1));
	    for ( int d = 0; d < target_dim; ++d )
	    {
		TEST_ASSERT( source_rank+1 == 
			     *(target_space_manager->field()->begin()
			       +n+d*num_points) );
	    }
	}
    }

    // Check the missing points.
    TEST_ASSERT( missing_points.size() > 0 );
    Teuchos::ArrayView<int> missed_in_map = 
	shared_domain_map.getMissedTargetPoints();
    TEST_ASSERT( missing_points.size() == missed_in_map.size() );

    std::sort( missing_points.begin(), missing_points.end() );
    std::sort( missed_in_map.begin(), missed_in_map.end() );

    for ( int n = 0; n < (int) missing_points.size(); ++n )
    {
	TEST_ASSERT( missing_points[n] == missed_in_map[n] );
    }
}

//---------------------------------------------------------------------------//
// Some points will be outside of the mesh in this test. The mesh is not
// rectilinear so we will be sure to search the kD-tree with points that
// aren't in the mesh.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_tiled_test10 )
{
    using namespace DataTransferKit;
    typedef MeshContainer<int> MeshType;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Compute element ordinal offsets so we make unique global ordinals.
    int edge_size = 10;
    int tri_offset_1 = 0;
    int quad_offset_1 = tri_offset_1 + (edge_size+1)*(edge_size+1)*2;
    int tri_offset_2 = quad_offset_1 + (edge_size+1)*(edge_size+1);
    int quad_offset_2 = tri_offset_2 + (edge_size+1)*(edge_size+1)*2;

    // Setup source mesh manager.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    if ( my_rank == 0 )
    {
	mesh_blocks[0] = 
	    buildTiledTriMesh( my_rank, my_size, edge_size, tri_offset_1 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else if ( my_rank == 1 )
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildTiledQuadMesh( my_rank, my_size, edge_size, quad_offset_1 );
    }
    else if ( my_rank == 2 )
    {
	mesh_blocks[0] = 
	    buildTiledTriMesh( my_rank, my_size, edge_size, tri_offset_2 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else 
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildTiledQuadMesh( my_rank, my_size, edge_size, quad_offset_2 );
    }
    comm->barrier();

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > source_mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 2 ) );

    // Setup target coordinate field manager.
    int point_dim = 2;
    Teuchos::RCP<MyField> coordinate_field =
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildTiledCoordinateField( my_rank, my_size, num_points, edge_size,
			       coordinate_field );
    Teuchos::RCP<FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<int,MyField> > source_evaluator;
    if ( my_rank == 0 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );
    }
    else if ( my_rank == 1 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[1], comm ) );
    }
    else if ( my_rank == 2 )
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );
    }
    else
    {
    	source_evaluator = Teuchos::rcp( new MyEvaluator( *mesh_blocks[1], comm ) );
    }
    comm->barrier();

    // Create data target. This target is a scalar.
    int target_dim = 1;
    Teuchos::RCP<MyField> target_field =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MeshType ,MyField> shared_domain_map( 
	comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data if it is in the mesh and 0.0 if it is outside.
    int source_rank;
    Teuchos::Array<int> missing_points;
    bool tagged;
    for ( int n = 0; n < num_points; ++n )
    {
	tagged = false;

	if ( *(coordinate_field->begin()+n) < 0.0 ||
	     *(coordinate_field->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field->begin()+n+num_points) > (edge_size-1)*my_size )
	{
	    missing_points.push_back(n);	
	    TEST_ASSERT( 0.0 == *(target_space_manager->field()->begin()+n) );
	    tagged = true;
	}
	
	else
	{
	    for ( int i = 0; i < my_size; ++i )
	    {
		if ( *(coordinate_field->begin()+n) >= (edge_size-1)*i &&
		     *(coordinate_field->begin()+n) <= (edge_size-1)*(i+1) &&
		     *(coordinate_field->begin()+n+num_points) >=
		     (edge_size-1)*i &&
		     *(coordinate_field->begin()+n+num_points) <= 
		     (edge_size-1)*(i+1) && !tagged )
		{
		    source_rank = std::floor(target_coord_manager->field()->getData()[n] 
					     / (edge_size-1));
		    TEST_ASSERT( source_rank+1 == 
				 target_space_manager->field()->getData()[n] );
		    tagged = true;
		}
	    }

	    if ( !tagged) 
	    {
		missing_points.push_back(n);	
		TEST_ASSERT( 0.0 == *(target_space_manager->field()->begin()+n) );
		tagged = true;
	    }
	}

	TEST_ASSERT( tagged );
    }

    // Check the missing points.
    TEST_ASSERT( missing_points.size() > 0 );
    Teuchos::ArrayView<int> missed_in_map = 
	shared_domain_map.getMissedTargetPoints();
    TEST_ASSERT( missing_points.size() == missed_in_map.size() );

    std::sort( missing_points.begin(), missing_points.end() );
    std::sort( missed_in_map.begin(), missed_in_map.end() );

    for ( int n = 0; n < (int) missing_points.size(); ++n )
    {
	TEST_ASSERT( missing_points[n] == missed_in_map[n] );
    }
}

//---------------------------------------------------------------------------//
// end tstSharedDomainMap10.cpp
//---------------------------------------------------------------------------//

