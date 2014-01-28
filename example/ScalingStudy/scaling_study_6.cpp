//---------------------------------------------------------------------------//
/*!
 * \file scaling_study_6.cpp
 * \author Stuart R. Slattery
 * \brief Scaling study 6. Strong scaling in 3D over cube of elements with
 * locally random target points.
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
// Mesh Implementation
//---------------------------------------------------------------------------//

class MyMesh
{
  public:

    typedef long int    global_ordinal_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<global_ordinal_type>& vertex_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<global_ordinal_type>& element_handles,
	    const Teuchos::Array<global_ordinal_type>& element_connectivity,
	    const Teuchos::Array<int>& permutation_list )
	: d_vertex_handles( vertex_handles )
	, d_coords( coords )
	, d_element_handles( element_handles )
	, d_element_connectivity( element_connectivity )
	, d_permutation_list( permutation_list )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    Teuchos::Array<global_ordinal_type>::const_iterator verticesBegin() const
    { return d_vertex_handles.begin(); }

    Teuchos::Array<global_ordinal_type>::const_iterator verticesEnd() const
    { return d_vertex_handles.end(); }

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
    
    Teuchos::Array<int>::const_iterator permutationBegin() const
    { return d_permutation_list.begin(); }

    Teuchos::Array<int>::const_iterator permutationEnd() const
    { return d_permutation_list.end(); }

  private:

    Teuchos::Array<global_ordinal_type> d_vertex_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<global_ordinal_type> d_element_handles;
    Teuchos::Array<global_ordinal_type> d_element_connectivity;
    Teuchos::Array<int> d_permutation_list;
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

    MyField( size_type size, int dim )
	: d_dim( dim )
	, d_data( size )
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
// Mesh traits specialization for MyMesh
template<>
class MeshTraits<MyMesh>
{
  public:

    typedef MyMesh::global_ordinal_type global_ordinal_type;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_vertex_iterator;
    typedef Teuchos::Array<double>::const_iterator 
    const_coordinate_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_element_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_connectivity_iterator;
    typedef Teuchos::Array<int>::const_iterator 
    const_permutation_iterator;

    static inline int vertexDim( const MyMesh& mesh )
    { return 3; }

    static inline const_vertex_iterator verticesBegin( const MyMesh& mesh )
    { return mesh.verticesBegin(); }

    static inline const_vertex_iterator verticesEnd( const MyMesh& mesh )
    { return mesh.verticesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline DTK_ElementTopology elementTopology( const MyMesh& mesh )
    { return DTK_HEXAHEDRON; }

    static inline int verticesPerElement( const MyMesh& mesh )
    { return 8; }


    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.elementsBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.elementsEnd(); }

    static inline const_connectivity_iterator 
    connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator 
    connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }

    static inline const_permutation_iterator 
    permutationBegin( const MyMesh& mesh )
    { return mesh.permutationBegin(); }

    static inline const_permutation_iterator 
    permutationEnd( const MyMesh& mesh )
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
class MyEvaluator : 
    public DataTransferKit::FieldEvaluator<MyMesh::global_ordinal_type,MyField>
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
Teuchos::RCP<MyMesh> buildMyMesh( int my_rank, int my_size, int edge_length )
{
    // Compute block indices.
    int num_blocks = std::pow( my_size, 1.0/3.0 );
    int k_block = std::floor( my_rank / (num_blocks*num_blocks) );
    int j_block = 
	std::floor( (my_rank-k_block*num_blocks*num_blocks) / num_blocks );
    int i_block = my_rank - j_block*num_blocks - k_block*num_blocks*num_blocks;

    // Make some vertices.
    int num_vertices = edge_length*edge_length*edge_length;
    int vertex_dim = 3;
    Teuchos::Array<long int> vertex_handles( num_vertices );
    Teuchos::Array<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int k = 0; k < edge_length; ++k )
    {
	for ( int j = 0; j < edge_length; ++j )
	{
	    for ( int i = 0; i < edge_length; ++i )
	    {
		idx = i + j*edge_length + k*edge_length*edge_length;
		vertex_handles[ idx ] = num_vertices*my_rank + idx;
		coords[ idx ] = i + i_block*(edge_length-1);
		coords[ num_vertices + idx ] = j + j_block*(edge_length-1);
		coords[ 2*num_vertices + idx ] = k + k_block*(edge_length-1);
	    }
	}
    }
    
    // Make the hexahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1)*(edge_length-1);
    Teuchos::Array<long int> hex_handles( num_elements );
    Teuchos::Array<long int> hex_connectivity( 8*num_elements );
    int elem_idx, vertex_idx;
    for ( int k = 0; k < (edge_length-1); ++k )
    {
	for ( int j = 0; j < (edge_length-1); ++j )
	{
	    for ( int i = 0; i < (edge_length-1); ++i )
	    {
		vertex_idx = i + j*edge_length + k*edge_length*edge_length;
		elem_idx = 
		    i + j*(edge_length-1) + k*(edge_length-1)*(edge_length-1);

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
		    = vertex_handles[vertex_idx+edge_length*edge_length];

		hex_connectivity[5*num_elements+elem_idx] 
		    = vertex_handles[vertex_idx+edge_length*edge_length+1];

		hex_connectivity[6*num_elements+elem_idx] 
		    = vertex_handles[
			vertex_idx+edge_length*edge_length+edge_length+1];

		hex_connectivity[7*num_elements+elem_idx] 
		    = vertex_handles[
			vertex_idx+edge_length*edge_length+edge_length];
	    }
	}
    }

    Teuchos::Array<int> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp(
	new MyMesh( vertex_handles, coords, hex_handles, hex_connectivity,
		    permutation_list ) );
}

//---------------------------------------------------------------------------//
// Coordinate field create function.
//---------------------------------------------------------------------------//
Teuchos::RCP<MyField> buildCoordinateField( int my_rank, int my_size, 
					    int num_points, int edge_size )
{
   // Compute block indices.
    int num_blocks = std::pow( my_size, 1.0/3.0 );
    int k_block = std::floor( my_rank / (num_blocks*num_blocks) );
    int j_block = 
	std::floor( (my_rank-k_block*num_blocks*num_blocks) / num_blocks );
    int i_block = my_rank - j_block*num_blocks - k_block*num_blocks*num_blocks;
    std::srand( my_rank*num_points*3 );
    int point_dim = 3;
    Teuchos::RCP<MyField> coordinate_field = 
	Teuchos::rcp( new MyField(num_points*point_dim, point_dim ) );

    for ( int i = 0; i < num_points; ++i )
    {
	*(coordinate_field->begin() + i) = 
	    (edge_size-1) * (double) std::rand() / RAND_MAX + 
	    (edge_size-1) * (num_blocks-i_block-1);
	*(coordinate_field->begin() + num_points + i ) = 
	    (edge_size-1) * (double) std::rand() / RAND_MAX + 
	    (edge_size-1) * (num_blocks-j_block-1);
	*(coordinate_field->begin() + 2*num_points + i ) = 
	    (edge_size-1) * (double) std::rand() / RAND_MAX + 
	    (edge_size-1) * (num_blocks-k_block-1);
    }

    return coordinate_field;
}

//---------------------------------------------------------------------------//
// Weak scaling study driver.
//---------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh.
    int global_size = 10;
    int edge_size = (global_size / std::pow(my_size,1.0/3.0) ) + 1;
    Teuchos::ArrayRCP<Teuchos::RCP<MyMesh> > mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP< MeshManager<MyMesh> > source_mesh_manager = Teuchos::rcp( 
	new MeshManager<MyMesh>( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field.
    int num_points = (edge_size-1)*(edge_size-1)*(edge_size-1);
    Teuchos::RCP<MyField> target_coords = 
	buildCoordinateField( my_rank, my_size, num_points, edge_size );
    Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	Teuchos::rcp( new FieldManager<MyField>( target_coords, comm ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<MyMesh::global_ordinal_type,MyField> > 
	source_evaluator = 
    	Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );

    // Create data target.
    int target_dim = 1;
    Teuchos::RCP<MyField> target_field = 
	Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP< FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );

    // Setup consistent interpolation mapping.
    SharedDomainMap<MyMesh,MyField> shared_domain_map( comm, 3 );

    // Setup the shared domain map ( this creates the mapping ).
    std::clock_t setup_start = clock();
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    std::clock_t setup_end = clock();

    // Apply the shared domain map ( this does the field evaluation and moves
    // the data ).
    std::clock_t apply_start = clock();
    shared_domain_map.apply( source_evaluator, target_space_manager );
    std::clock_t apply_end = clock();

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data. Count the number of target points that
    // failed the test.
    comm->barrier();
    int num_blocks = std::pow( my_size, 1.0/3.0 );
    int k_block = std::floor( my_rank / (num_blocks*num_blocks) );
    int j_block = 
	std::floor( (my_rank-k_block*num_blocks*num_blocks) / num_blocks );
    int i_block = my_rank - j_block*num_blocks - k_block*num_blocks*num_blocks;
    int source_rank = (num_blocks-i_block-1) + (num_blocks-j_block-1)*num_blocks +
		      (num_blocks-k_block-1)*num_blocks*num_blocks;
    int local_test_failed = 0;
    for ( long int n = 0; n < target_space_manager->field()->size(); ++n )
    {
    	if ( source_rank+1 != target_space_manager->field()->getData()[n] )
    	{
    	    local_test_failed += 1;
    	}
    }
    comm->barrier();

    int global_test_failed;
    Teuchos::reduceAll<int,int>( *comm,
    				 Teuchos::REDUCE_SUM,
    				 1,
    				 &local_test_failed,
    				 &global_test_failed );
   
    if ( my_rank == 0 )
    {
    	if ( global_test_failed )
    	{
    	    std::cout << std::endl << "TEST FAILURE " 
		      << global_test_failed << std::endl;
    	}
    	else
    	{
    	    std::cout << std::endl << "TEST PASSED" << std::endl;
    	}
    }
    comm->barrier();

    // Timing.
    double local_setup_time = 
    	(double)(setup_end - setup_start) / CLOCKS_PER_SEC;

    double global_min_setup_time;
    Teuchos::reduceAll<int,double>( *comm,
    				    Teuchos::REDUCE_MIN,
    				    1,
    				    &local_setup_time,
    				    &global_min_setup_time );

    double global_max_setup_time;
    Teuchos::reduceAll<int,double>( *comm,
    				    Teuchos::REDUCE_MAX,
    				    1,
    				    &local_setup_time,
    				    &global_max_setup_time );

    double global_average_setup_time;
    Teuchos::reduceAll<int,double>( *comm,
    				    Teuchos::REDUCE_SUM,
    				    1,
    				    &local_setup_time,
    				    &global_average_setup_time );
    global_average_setup_time /= my_size;

    double local_apply_time = 
    	(double)(apply_end - apply_start) / CLOCKS_PER_SEC;

    double global_min_apply_time;
    Teuchos::reduceAll<int,double>( *comm,
    				    Teuchos::REDUCE_MIN,
    				    1,
    				    &local_apply_time,
    				    &global_min_apply_time );

    double global_max_apply_time;
    Teuchos::reduceAll<int,double>( *comm,
    				    Teuchos::REDUCE_MAX,
    				    1,
    				    &local_apply_time,
    				    &global_max_apply_time );

    double global_average_apply_time;
    Teuchos::reduceAll<int,double>( *comm,
    				    Teuchos::REDUCE_SUM,
    				    1,
    				    &local_apply_time,
    				    &global_average_apply_time );
    global_average_apply_time /= my_size;

    comm->barrier();

    if ( my_rank == 0 )
    {
    	std::cout << "==================================================" 
    		  << std::endl;
    	std::cout << "DTK weak scaling study" << std::endl;
    	std::cout << "Number of processors:      " << my_size << std::endl;
    	std::cout << "Local number of elements:  " 
		  << (edge_size-1)*(edge_size-1)*(edge_size-1)
    		  << std::endl;
    	std::cout << "Local number of points:    " 
		  << (edge_size-1)*(edge_size-1)*(edge_size-1)
    		  << std::endl;
    	std::cout << "Global number of elements: " 
		  << (edge_size-1)*(edge_size-1)*(edge_size-1)*my_size 
    		  << std::endl;
    	std::cout << "Global number of points:   " 
		  << (edge_size-1)*(edge_size-1)*(edge_size-1)*my_size 
    		  << std::endl;
    	std::cout << "--------------------------------------------------"
    		  << std::endl;
    	std::cout << "Global min setup time (s):     " 
    		  << global_min_setup_time << std::endl;
    	std::cout << "Global max setup time (s):     " 
    		  << global_max_setup_time << std::endl;
    	std::cout << "Global average setup time (s): " 
    		  << global_average_setup_time << std::endl;
    	std::cout << "--------------------------------------------------"
    		  << std::endl;
    	std::cout << "Global min apply time (s):     " 
    		  << global_min_apply_time << std::endl;
    	std::cout << "Global max apply time (s):     " 
    		  << global_max_apply_time << std::endl;
    	std::cout << "Global average apply time (s): " 
    		  << global_average_apply_time << std::endl;
    	std::cout << "==================================================" 
    		  << std::endl;
    }

    comm->barrier();

    return 0;
}

//---------------------------------------------------------------------------//
// end scaling_study_6.cpp
//---------------------------------------------------------------------------//

