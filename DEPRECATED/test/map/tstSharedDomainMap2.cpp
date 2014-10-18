//---------------------------------------------------------------------------//
/*!
 * \file tstSharedDomainMap2.cpp
 * \author Stuart R. Slattery
 * \brief Shared domain map unit test 2.
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
#include <DTK_MeshBlock.hpp>
#include <DTK_MeshManager.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
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
	, d_data( dim*size )
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

    Teuchos::ArrayView<double> getData()
    { return d_data(); }

    const Teuchos::ArrayView<const double> getData() const
    { return d_data(); }

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
// FieldEvaluator Implementation. This one populates a 3-vector.
class MyEvaluator : 
    public DataTransferKit::FieldEvaluator<DataTransferKit::MeshId,MyField>
{
  public:

    MyEvaluator( const Teuchos::RCP<DataTransferKit::MeshBlock>& mesh, 
		 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
	: d_mesh( mesh )
	, d_comm( comm )
    { /* ... */ }

    ~MyEvaluator()
    { /* ... */ }

    MyField evaluate( 
	const Teuchos::ArrayRCP<DataTransferKit::MeshId>& elements,
	const Teuchos::ArrayRCP<double>& coords )
    {
	int num_elements = elements.size();
	MyField evaluated_data( num_elements, 3 );
	for ( int n = 0; n < num_elements; ++n )
	{
	    if ( std::find( d_mesh->elementIds().begin(),
			    d_mesh->elementIds().end(),
			    elements[n] ) != d_mesh->elementIds().end() )
	    {
		*(evaluated_data.begin() + n ) = d_comm->getRank() + 1.0;
		*(evaluated_data.begin() + num_elements + n ) = 
		    d_comm->getRank() + 1.0;
		*(evaluated_data.begin() + 2*num_elements + n ) = 
		    d_comm->getRank() + 1.0;
	    }
	    else
	    {
 		*(evaluated_data.begin() + n ) = 0.0;
		*(evaluated_data.begin() + num_elements + n ) = 0.0;
		*(evaluated_data.begin() + 2*num_elements + n ) = 0.0;
	    }
	}
	return evaluated_data;
    }

  private:

    Teuchos::RCP<DataTransferKit::MeshBlock> d_mesh;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// Mesh create functions.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshBlock> 
buildMyMesh( int my_rank, int my_size, int edge_length )
{
    // Make some vertices.
    int num_vertices = edge_length*edge_length*2;
    int vertex_dim = 3;
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handles( num_vertices );
    Teuchos::ArrayRCP<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    vertex_handles[ idx ] = (DataTransferKit::MeshId) num_vertices*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j;
	    coords[ 2*num_vertices + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + num_vertices / 2;
	    vertex_handles[ idx ] = (DataTransferKit::MeshId) num_vertices*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j;
	    coords[ 2*num_vertices + idx ] = 1.0;
	}
    }
    
    // Make the hexahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::ArrayRCP<DataTransferKit::MeshId> hex_handles( num_elements );
    Teuchos::ArrayRCP<DataTransferKit::MeshId> hex_connectivity( 8*num_elements );
    int elem_idx, vertex_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    vertex_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    hex_handles[elem_idx] = num_elements*my_rank + elem_idx;

	    hex_connectivity[elem_idx] 
		= vertex_handles[vertex_idx];

	    hex_connectivity[num_elements+elem_idx] 
		= vertex_handles[vertex_idx+1];

	    hex_connectivity[2*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length+1];

	    hex_connectivity[3*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length];

	    hex_connectivity[4*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2];

	    hex_connectivity[5*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2+1];

 	    hex_connectivity[6*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2+edge_length+1];

	    hex_connectivity[7*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2+edge_length];
	}
    }

    Teuchos::ArrayRCP<int> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp( 
	new DataTransferKit::MeshContainer( 
	    3, vertex_handles, coords, 
	    DataTransferKit::DTK_HEXAHEDRON, 8, 
	    hex_handles, hex_connectivity,
	    permutation_list ) );
}

//---------------------------------------------------------------------------//
// This function creates a tiled mesh such that the global mesh is not
// rectilinear.
Teuchos::RCP<DataTransferKit::MeshBlock> 
buildTiledMesh( int my_rank, int my_size, int edge_length )
{
    // Make some vertices.
    int num_vertices = edge_length*edge_length*2;
    int vertex_dim = 3;
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handles( num_vertices );
    Teuchos::ArrayRCP<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    vertex_handles[ idx ] = (DataTransferKit::MeshId) num_vertices*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j + my_rank*(edge_length-1);
	    coords[ 2*num_vertices + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + num_vertices / 2;
	    vertex_handles[ idx ] = (DataTransferKit::MeshId) num_vertices*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_vertices + idx ] = j + my_rank*(edge_length-1);
	    coords[ 2*num_vertices + idx ] = 1.0;
	}
    }
    
    // Make the hexahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::ArrayRCP<DataTransferKit::MeshId> hex_handles( num_elements );
    Teuchos::ArrayRCP<DataTransferKit::MeshId> hex_connectivity( 8*num_elements );
    int elem_idx, vertex_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    vertex_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    hex_handles[elem_idx] = num_elements*my_rank + elem_idx;

	    hex_connectivity[elem_idx] 
		= vertex_handles[vertex_idx];

	    hex_connectivity[num_elements+elem_idx] 
		= vertex_handles[vertex_idx+1];

	    hex_connectivity[2*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length+1];

	    hex_connectivity[3*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+edge_length];

	    hex_connectivity[4*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2];

	    hex_connectivity[5*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2+1];

 	    hex_connectivity[6*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2+edge_length+1];

	    hex_connectivity[7*num_elements+elem_idx] 
		= vertex_handles[vertex_idx+num_vertices/2+edge_length];
	}
    }

    Teuchos::ArrayRCP<int> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp(
	new DataTransferKit::MeshContainer( 
	    3,vertex_handles, coords, 
	    DataTransferKit::DTK_HEXAHEDRON, 8, 
	    hex_handles, hex_connectivity,
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
	*(coordinate_field->begin() + 2*num_points + i ) = 
	    (double) std::rand() / RAND_MAX;
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
	*(coordinate_field->begin() + 2*num_points + i ) = 
	    1.2 * (double) std::rand() / RAND_MAX - 0.1;
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
	*(coordinate_field->begin() + 2*num_points + i ) = 
	    1.2 * (double) std::rand() / RAND_MAX - 0.1;
    }
}

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
// All points will be in the mesh in this test.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_test2 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP< MeshManager > source_mesh_manager = Teuchos::rcp( 
	new MeshManager( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field manager.
    int num_points = 1000;
    int point_dim = 3;
    Teuchos::RCP<MyField> coordinate_field =
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildCoordinateField( my_rank, my_size, num_points, edge_size,
			  coordinate_field );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<MeshId,MyField> > source_evaluator = 
    	Teuchos::rcp( new MyEvaluator( mesh_blocks[0], comm ) );

    // Create data target. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );
    // Setup and apply the shared domain mapping.
    SharedDomainMap<MyField> shared_domain_map( 
	comm, source_mesh_manager->dim() );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data.
    int source_rank;
    for ( int n = 0; n < num_points; ++n )
    {
	source_rank = std::floor(*(coordinate_field->begin()+n) / (edge_size-1));
	TEST_EQUALITY( source_rank+1,
		       *(target_space_manager->field()->begin()+n) );
	TEST_EQUALITY( source_rank+1, 
		       *(target_space_manager->field()->begin()+n+num_points) );
	TEST_EQUALITY( source_rank+1, 
		       *(target_space_manager->field()->begin()+n+2*num_points) );
    }
}

//---------------------------------------------------------------------------//
// Some points will be outside of the mesh in this test.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_expanded_test2 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP< MeshManager > source_mesh_manager = Teuchos::rcp( 
	new MeshManager( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field manager.
    int num_points = 1000;
    int point_dim = 3;
    Teuchos::RCP<MyField> coordinate_field =
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
				  coordinate_field );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<MeshId,MyField> > 
	source_evaluator = 
    	Teuchos::rcp( new MyEvaluator( mesh_blocks[0], comm ) );

    // Create data target. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MyField> shared_domain_map( 
	comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data if it is in the mesh and 0.0 if it is outside.
    int source_rank;
    Teuchos::Array<MeshId> missing_points;
    for ( int n = 0; n < num_points; ++n )
    {
	if ( *(coordinate_field->begin()+n) < 0.0 ||
	     *(coordinate_field->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field->begin()+n+num_points) > edge_size-1 ||
	     *(coordinate_field->begin()+n+2*num_points) < 0.0 ||
	     *(coordinate_field->begin()+n+2*num_points) > 1.0 )
	{
	    missing_points.push_back(n);	
	    TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()+n) );
	    TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()
				  +n+num_points) );
	    TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()
				  +n+2*num_points) );
	}
	else
	{
	    source_rank = std::floor(target_coord_manager->field()->getData()[n] 
	    			     / (edge_size-1));
	    TEST_EQUALITY( source_rank+1, 
			   target_space_manager->field()->getData()[n] );
	    TEST_EQUALITY( source_rank+1, 
			   target_space_manager->field()->getData()[
			       n + num_points] );
	    TEST_EQUALITY( source_rank+1,
			   target_space_manager->field()->getData()[
			       n + 2*num_points] );
	}
    }

    // Check the missing points.
    TEST_ASSERT( missing_points.size() > 0 );
    Teuchos::ArrayView<MeshId> missed_in_map = 
	shared_domain_map.getMissedTargetPoints();
    TEST_EQUALITY( missing_points.size(), missed_in_map.size() );

    std::sort( missing_points.begin(), missing_points.end() );
    std::sort( missed_in_map.begin(), missed_in_map.end() );

    for ( int n = 0; n < (int) missing_points.size(); ++n )
    {
	TEST_EQUALITY( missing_points[n], missed_in_map[n] );
    }
}

//---------------------------------------------------------------------------//
// Some points will be outside of the mesh in this test. The mesh is not
// rectilinear so we will be sure to search the kD-tree with points that
// aren't in the mesh.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_tiled_test2 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTiledMesh( my_rank, my_size, edge_size );
    Teuchos::RCP< MeshManager > source_mesh_manager = Teuchos::rcp( 
	new MeshManager( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field manager.
    int num_points = 1000;
    int point_dim = 3;
    Teuchos::RCP<MyField> coordinate_field =
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildTiledCoordinateField( my_rank, my_size, num_points, edge_size,
			       coordinate_field );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<MeshId,MyField> > source_evaluator = 
    	Teuchos::rcp( new MyEvaluator( mesh_blocks[0], comm ) );

    // Create data target. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MyField> shared_domain_map( 
	comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data if it is in the mesh and 0.0 if it is outside.
    int source_rank;
    Teuchos::Array<MeshId> missing_points;
    bool tagged;
    for ( int n = 0; n < num_points; ++n )
    {
	tagged = false;

	if ( *(coordinate_field->begin()+n) < 0.0 ||
	     *(coordinate_field->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field->begin()+n+num_points) > (edge_size-1)*my_size ||
	     *(coordinate_field->begin()+n+2*num_points) < 0.0 ||
	     *(coordinate_field->begin()+n+2*num_points) > 1.0 )
	{
	    missing_points.push_back(n);	
	    TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()+n) );
	    TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()
				  +n+num_points) );
	    TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()
				  +n+2*num_points) );
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
		     (edge_size-1)*(i+1) &&
		     *(coordinate_field->begin()+n+2*num_points) >= 0.0 &&
		     *(coordinate_field->begin()+n+2*num_points) <= 1.0 && !tagged )
		{
		    source_rank = std::floor(target_coord_manager->field()->getData()[n] 
					     / (edge_size-1));
		    TEST_EQUALITY( source_rank+1,
				   target_space_manager->field()->getData()[n] );
		    TEST_EQUALITY( source_rank+1,
				   target_space_manager->field()->getData()[
				       n + num_points] );
		    TEST_EQUALITY( source_rank+1,
				   target_space_manager->field()->getData()[
				       n + 2*num_points] );
		    tagged = true;
		}
	    }

	    if ( !tagged) 
	    {
		missing_points.push_back(n);	
		TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()+n) );
		TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()
				      +n+num_points) );
		TEST_EQUALITY( 0.0, *(target_space_manager->field()->begin()
				      +n+2*num_points) );
		tagged = true;
	    }
	}

	TEST_ASSERT( tagged );
    }

    // Check the missing points.
    TEST_ASSERT( missing_points.size() > 0 );
    Teuchos::ArrayView<MeshId> missed_in_map = 
	shared_domain_map.getMissedTargetPoints();
    TEST_EQUALITY( missing_points.size(), missed_in_map.size() );

    std::sort( missing_points.begin(), missing_points.end() );
    std::sort( missed_in_map.begin(), missed_in_map.end() );

    for ( int n = 0; n < (int) missing_points.size(); ++n )
    {
	TEST_EQUALITY( missing_points[n], missed_in_map[n] );
    }
}

//---------------------------------------------------------------------------//
// This test has one source mesh and one field evaluator for that source mesh
// but two target coordinates and target spaces that it will be mapped to. We
// do this to make sure that we can reuse a mesh description with multiple
// targets.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_two_maps_test2 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP< MeshManager > source_mesh_manager = Teuchos::rcp( 
	new MeshManager( mesh_blocks, comm, 3 ) );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<MeshId,MyField> >
	source_evaluator = 
    	Teuchos::rcp( new MyEvaluator( mesh_blocks[0], comm ) );

    int num_points = 1000;
    int point_dim = 3;

    // Setup target coordinate field manager 1.
    Teuchos::RCP<MyField> coordinate_field_1 =
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
				  coordinate_field_1 );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager_1 = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field_1, comm ) );

    // Create data target 1. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field_1 =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager_1 = 
	Teuchos::rcp( 
	new FieldManager<MyField>( target_field_1, comm ) );

    // Setup target coordinate field manager 2.
    Teuchos::RCP<MyField> coordinate_field_2 = 
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
				  coordinate_field_2 );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager_2 = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field_2, comm ) );


    // Create data target 2. This target is a 3-vector.
    Teuchos::RCP<MyField> target_field_2 =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager_2 = 
	Teuchos::rcp( 
	new FieldManager<MyField>( target_field_2, comm ) );

    // Setup and apply the shared domain mapping 1.
    SharedDomainMap<MyField> shared_domain_map_1( 
	comm, source_mesh_manager->dim(), true );
    shared_domain_map_1.setup( source_mesh_manager, target_coord_manager_1 );
    shared_domain_map_1.apply( source_evaluator, target_space_manager_1 );

    // Check the first data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.  
    int source_rank;
    Teuchos::Array<MeshId> missing_points_1;
    for ( int n = 0; n < num_points; ++n )
    {
	if ( *(coordinate_field_1->begin()+n) < 0.0 ||
	     *(coordinate_field_1->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field_1->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field_1->begin()+n+num_points) > edge_size-1 ||
	     *(coordinate_field_1->begin()+n+2*num_points) < 0.0 ||
	     *(coordinate_field_1->begin()+n+2*num_points) > 1.0 )
	{
	    missing_points_1.push_back(n);	
	    TEST_ASSERT( 0.0 == *(target_space_manager_1->field()->begin()+n) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_1->field()->begin()
				  +n+num_points) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_1->field()->begin()
				  +n+2*num_points) );
	}
	else
	{
	    source_rank = std::floor(
		target_coord_manager_1->field()->getData()[n] / (edge_size-1));
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_1->field()->getData()[n] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_1->field()->getData()[
	    		     n + num_points] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_1->field()->getData()[
	    		     n + 2*num_points] );
	}
    }

    // Check the missing points.
    TEST_ASSERT( missing_points_1.size() > 0 );
    Teuchos::ArrayView<MeshId> missed_in_map_1 = 
	shared_domain_map_1.getMissedTargetPoints();
    TEST_ASSERT( missing_points_1.size() == missed_in_map_1.size() );

    std::sort( missing_points_1.begin(), missing_points_1.end() );
    std::sort( missed_in_map_1.begin(), missed_in_map_1.end() );

    for ( int n = 0; n < (int) missing_points_1.size(); ++n )
    {
	TEST_ASSERT( missing_points_1[n] == missed_in_map_1[n] );
    }

    // Setup and apply the shared domain mapping 2.
    SharedDomainMap<MyField> shared_domain_map_2( 
	comm, source_mesh_manager->dim(), true );
    shared_domain_map_2.setup( source_mesh_manager, target_coord_manager_2 );
    shared_domain_map_2.apply( source_evaluator, target_space_manager_2 );

    // Check the second data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.  int source_rank;
    Teuchos::Array<MeshId> missing_points_2;
    for ( int n = 0; n < num_points; ++n )
    {
	if ( *(coordinate_field_2->begin()+n) < 0.0 ||
	     *(coordinate_field_2->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field_2->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field_2->begin()+n+num_points) > edge_size-1 ||
	     *(coordinate_field_2->begin()+n+2*num_points) < 0.0 ||
	     *(coordinate_field_2->begin()+n+2*num_points) > 1.0 )
	{
	    missing_points_2.push_back(n);	
	    TEST_ASSERT( 0.0 == *(target_space_manager_2->field()->begin()+n) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_2->field()->begin()
				  +n+num_points) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_2->field()->begin()
				  +n+2*num_points) );
	}
	else
	{
	    source_rank = std::floor(
		target_coord_manager_2->field()->getData()[n] / (edge_size-1));
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_2->field()->getData()[n] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_2->field()->getData()[
	    		     n + num_points] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_2->field()->getData()[
	    		     n + 2*num_points] );
	}
    }

    // Check the missing points 2.
    TEST_ASSERT( missing_points_2.size() > 0 );
    Teuchos::ArrayView<MeshId> missed_in_map_2 = 
	shared_domain_map_2.getMissedTargetPoints();
    TEST_ASSERT( missing_points_2.size() == missed_in_map_2.size() );

    std::sort( missing_points_2.begin(), missing_points_2.end() );
    std::sort( missed_in_map_2.begin(), missed_in_map_2.end() );

    for ( int n = 0; n < (int) missing_points_2.size(); ++n )
    {
	TEST_ASSERT( missing_points_2[n] == missed_in_map_2[n] );
    }
}

//---------------------------------------------------------------------------//
// This test has one source mesh and one field evaluator for that source mesh
// but two target coordinates and target spaces that it will be mapped to. We
// do this to make sure that we can use multiple targets with a single map.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_two_targets_test2 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP<MeshManager > source_mesh_manager = Teuchos::rcp( 
	new MeshManager( mesh_blocks, comm, 3 ) );

    // Create field evaluator.
    Teuchos::RCP< FieldEvaluator<MeshId,MyField> > 
	source_evaluator = 
    	Teuchos::rcp( new MyEvaluator( mesh_blocks[0], comm ) );

    int num_points = 1000;
    int point_dim = 3;

    // Setup target coordinate field manager 1.
    Teuchos::RCP<MyField> coordinate_field_1 = 
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
				  coordinate_field_1 );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager_1 = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field_1, comm ) );

    // Create data target 1. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field_1 =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager_1 = 
	Teuchos::rcp( 
	new FieldManager<MyField>( target_field_1, comm ) );

    // Setup target coordinate field manager 2.
    Teuchos::RCP<MyField> coordinate_field_2 =
	Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
				  coordinate_field_2 );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager_2 = 
	Teuchos::rcp( new FieldManager<MyField>( coordinate_field_2, comm ) );


    // Create data target 2. This target is a 3-vector.
    Teuchos::RCP<MyField> target_field_2 =  
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager_2 = 
	Teuchos::rcp( 
	new FieldManager<MyField>( target_field_2, comm ) );

    // Setup shared domain mapping.
    SharedDomainMap<MyField> shared_domain_map( 
	comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager_1 );

    // Apply the shared domain mapping 1.
    shared_domain_map.apply( source_evaluator, target_space_manager_1 );

    // Check the first data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.  
    int source_rank;
    Teuchos::Array<MeshId> missing_points_1;
    for ( int n = 0; n < num_points; ++n )
    {
	if ( *(coordinate_field_1->begin()+n) < 0.0 ||
	     *(coordinate_field_1->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field_1->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field_1->begin()+n+num_points) > edge_size-1 ||
	     *(coordinate_field_1->begin()+n+2*num_points) < 0.0 ||
	     *(coordinate_field_1->begin()+n+2*num_points) > 1.0 )
	{
	    missing_points_1.push_back(n);	
	    TEST_ASSERT( 0.0 == *(target_space_manager_1->field()->begin()+n) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_1->field()->begin()
				  +n+num_points) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_1->field()->begin()
				  +n+2*num_points) );
	}
	else
	{
	    source_rank = std::floor(
		target_coord_manager_1->field()->getData()[n] / (edge_size-1));
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_1->field()->getData()[n] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_1->field()->getData()[
	    		     n + num_points] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_1->field()->getData()[
	    		     n + 2*num_points] );
	}
    }

    // Check the missing points.
    TEST_ASSERT( missing_points_1.size() > 0 );
    Teuchos::ArrayView<MeshId> missed_in_map_1 = 
	shared_domain_map.getMissedTargetPoints();
    TEST_ASSERT( missing_points_1.size() == missed_in_map_1.size() );

    std::sort( missing_points_1.begin(), missing_points_1.end() );
    std::sort( missed_in_map_1.begin(), missed_in_map_1.end() );

    for ( int n = 0; n < (int) missing_points_1.size(); ++n )
    {
	TEST_ASSERT( missing_points_1[n] == missed_in_map_1[n] );
    }

    // Apply the shared domain mapping 2.
    shared_domain_map.apply( source_evaluator, target_space_manager_2 );

    // Check the second data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.  int source_rank;
    Teuchos::Array<MeshId> missing_points_2;
    for ( int n = 0; n < num_points; ++n )
    {
	if ( *(coordinate_field_2->begin()+n) < 0.0 ||
	     *(coordinate_field_2->begin()+n) > (edge_size-1)*my_size ||
	     *(coordinate_field_2->begin()+n+num_points) < 0.0 ||
	     *(coordinate_field_2->begin()+n+num_points) > edge_size-1 ||
	     *(coordinate_field_2->begin()+n+2*num_points) < 0.0 ||
	     *(coordinate_field_2->begin()+n+2*num_points) > 1.0 )
	{
	    missing_points_2.push_back(n);	
	    TEST_ASSERT( 0.0 == *(target_space_manager_2->field()->begin()+n) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_2->field()->begin()
				  +n+num_points) );
	    TEST_ASSERT( 0.0 == *(target_space_manager_2->field()->begin()
				  +n+2*num_points) );
	}
	else
	{
	    source_rank = std::floor(
		target_coord_manager_2->field()->getData()[n] / (edge_size-1));
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_2->field()->getData()[n] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_2->field()->getData()[
	    		     n + num_points] );
	    TEST_ASSERT( source_rank+1 == 
	    		 target_space_manager_2->field()->getData()[
	    		     n + 2*num_points] );
	}
    }

    // Check the missing points 2.
    TEST_ASSERT( missing_points_2.size() > 0 );
    Teuchos::ArrayView<MeshId> missed_in_map_2 = 
	shared_domain_map.getMissedTargetPoints();
    TEST_ASSERT( missing_points_2.size() == missed_in_map_2.size() );

    std::sort( missing_points_2.begin(), missing_points_2.end() );
    std::sort( missed_in_map_2.begin(), missed_in_map_2.end() );

    for ( int n = 0; n < (int) missing_points_2.size(); ++n )
    {
	TEST_ASSERT( missing_points_2[n] == missed_in_map_2[n] );
    }
}

//---------------------------------------------------------------------------//
// end tstSharedDomainMap2.cpp
//---------------------------------------------------------------------------//

