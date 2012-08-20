//---------------------------------------------------------------------------//
/*!
 * \file tstSharedDomainMap1.cpp
 * \author Stuart R. Slattery
 * \brief Shared domain map unit tests 1.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_SharedDomainMap.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshManager.hpp>

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

    MyMesh( const Teuchos::Array<int>& vertex_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<int>& quad_handles,
	    const Teuchos::Array<int>& quad_connectivity,
	    const Teuchos::Array<int>& permutation_list )
	: d_vertex_handles( vertex_handles )
	, d_coords( coords )
	, d_quad_handles( quad_handles )
	, d_quad_connectivity( quad_connectivity )
	, d_permutation_list( permutation_list )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    Teuchos::Array<int>::const_iterator verticesBegin() const
    { return d_vertex_handles.begin(); }

    Teuchos::Array<int>::const_iterator verticesEnd() const
    { return d_vertex_handles.end(); }

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
    
    Teuchos::Array<int>::const_iterator permutationBegin() const
    { return d_permutation_list.begin(); }

    Teuchos::Array<int>::const_iterator permutationEnd() const
    { return d_permutation_list.end(); }


  private:

    Teuchos::Array<int> d_vertex_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<int> d_quad_handles;
    Teuchos::Array<int> d_quad_connectivity;
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
    typedef Teuchos::Array<int>::const_iterator const_vertex_iterator;
    typedef Teuchos::Array<double>::const_iterator const_coordinate_iterator;
    typedef Teuchos::Array<int>::const_iterator const_element_iterator;
    typedef Teuchos::Array<int>::const_iterator const_connectivity_iterator;
    typedef Teuchos::Array<int>::const_iterator 
    const_permutation_iterator;

    static inline int vertexDim( const MyMesh& mesh )
    { return 2; }

    static inline const_vertex_iterator verticesBegin( const MyMesh& mesh )
    { return mesh.verticesBegin(); }

    static inline const_vertex_iterator verticesEnd( const MyMesh& mesh )
    { return mesh.verticesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline DTK_ElementTopology elementTopology( const MyMesh& mesh )
    { return DTK_QUADRILATERAL; }

    static inline int verticesPerElement( const MyMesh& mesh )
    { return 4; }


    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.quadsBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.quadsEnd(); }

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
// Mesh create funciton.
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
Teuchos::RCP<MyMesh> buildMyMesh()
{
    int my_rank = getDefaultComm<int>()->getRank();

    // Make some vertices.
    int num_vertices = 10;
    int vertex_dim = 2;
    Teuchos::Array<int> vertex_handles( num_vertices );
    Teuchos::Array<double> coords( vertex_dim*num_vertices );

    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles[i] = (num_vertices / 2)*my_rank + i;
    }
    for ( int i = 0; i < num_vertices / 2; ++i )
    {
	coords[ i ] = my_rank;
	coords[ num_vertices + i ] = i;
    }
    for ( int i = num_vertices / 2; i < num_vertices; ++i )
    {
	coords[ i ] = my_rank + 1;
	coords[ num_vertices + i ] = i - num_vertices/2;
    }
    
    // Make the quads.
    int num_quads = 4;
    Teuchos::Array<int> quad_handles( num_quads );
    Teuchos::Array<int> quad_connectivity( 4*num_quads );
    
    for ( int i = 0; i < num_quads; ++i )
    {
	quad_handles[ i ] = num_quads*my_rank + i;
	quad_connectivity[ i ] = vertex_handles[i];
	quad_connectivity[ num_quads + i ] = vertex_handles[num_vertices/2 + i];
	quad_connectivity[ 2*num_quads + i ] = vertex_handles[num_vertices/2 + i + 1];
	quad_connectivity[ 3*num_quads + i ] = vertex_handles[i+1];
    }

    Teuchos::Array<int> permutation_list( 4 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return Teuchos::rcp( 
	new MyMesh( vertex_handles, coords, quad_handles, quad_connectivity,
		    permutation_list ) );
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
Teuchos::RCP<MyField> buildCoordinateField()
{
    int num_points = 4;
    int point_dim = 2;
    Teuchos::RCP<MyField> coordinate_field = 
	Teuchos::rcp( new MyField( num_points*point_dim, point_dim ) );

    for ( int i = 0; i < num_points; ++i )
    {
	*(coordinate_field->begin() + i) = i + 0.5;
	*(coordinate_field->begin() + num_points + i ) = 
	    getDefaultComm<int>()->getRank() + 0.5;
    }

    return coordinate_field;
}

//---------------------------------------------------------------------------//
// Unit tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_size = comm->getSize();

    // This is a 4 processor test.
    if ( my_size == 4 )
    {
	// Setup source mesh manager.
	Teuchos::ArrayRCP<Teuchos::RCP<MyMesh> > mesh_blocks( 1 );
	mesh_blocks[0] = buildMyMesh();
	Teuchos::RCP< MeshManager<MyMesh> > source_mesh_manager = Teuchos::rcp(
	    new MeshManager<MyMesh>( mesh_blocks, comm, 2 ) );

	// Setup target coordinate field manager
	Teuchos::RCP< FieldManager<MyField> > target_coord_manager = 
	    Teuchos::rcp( 
		new FieldManager<MyField>( buildCoordinateField(), comm ) );

	// Create field evaluator.
	Teuchos::RCP< FieldEvaluator<MyMesh,MyField> > source_evaluator = 
	    Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );

	// Create data target manager
	int field_size = target_coord_manager->field()->size() 
			 / target_coord_manager->field()->dim();
	Teuchos::RCP<MyField> target_field = 
	    Teuchos::rcp( new MyField( field_size, 1 ) );
	Teuchos::RCP< FieldManager<MyField> > target_space_manager = 
	    Teuchos::rcp( 
		new FieldManager<MyField>( target_field, comm ) );

	// Setup and apply the evaluation to the field.
	SharedDomainMap<MyMesh,MyField> shared_domain_map( 
	    comm, source_mesh_manager->dim() );
	shared_domain_map.setup( source_mesh_manager, target_coord_manager );
	shared_domain_map.apply( source_evaluator, target_space_manager );

	// Check the data transfer.
	for ( int n = 0; n < target_space_manager->field()->size(); ++n )
	{
	    TEST_ASSERT( *(target_space_manager->field()->begin()+n) 
			 == n + 1 );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstSharedDomainMap1.cpp
//---------------------------------------------------------------------------//
