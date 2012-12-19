//---------------------------------------------------------------------------//
/*!
 * \file tstIntegralAssemblyMap6.cpp
 * \author Stuart R. Slattery
 * \brief Integral assembly map unit test 6 for repeated geometries.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include <DTK_IntegralAssemblyMap.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldIntegrator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_GeometryTraits.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Ptr.hpp>
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
// FieldIntegrator Implementation.
class MyIntegrator : public DataTransferKit::FieldIntegrator<
    DataTransferKit::MeshContainer<int> ,MyField>
{
  public:

    MyIntegrator( const DataTransferKit::MeshContainer<int>& mesh, 
		 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
	: d_mesh( mesh )
	, d_comm( comm )
    { /* ... */ }

    ~MyIntegrator()
    { /* ... */ }

    // If the global id is valid, then set the element integral to 2.0
    MyField integrate( 
	const Teuchos::ArrayRCP<
	    DataTransferKit::MeshContainer<int>::global_ordinal_type>& elements )
    {
	int num_elements = elements.size();
	MyField integrated_data( num_elements, 3 );
	for ( int n = 0; n < elements.size(); ++n )
	{
	    if ( std::find( d_mesh.elementsBegin(),
			    d_mesh.elementsEnd(),
			    elements[n] ) != d_mesh.elementsEnd() )
	    {
		*(integrated_data.begin() + n ) = 2.0;
		*(integrated_data.begin() + n + num_elements) = 2.0;
		*(integrated_data.begin() + n + 2*num_elements) = 2.0;
	    }
	    else
	    {
 		*(integrated_data.begin() + n ) = 6789.443;
		*(integrated_data.begin() + n + num_elements) = 6789.443;
		*(integrated_data.begin() + n + 2*num_elements) = 6789.443;
	    }
	}
	return integrated_data;
    }

  private:

    DataTransferKit::MeshContainer<int>  d_mesh;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// ElementMeasure Implementation.
class MyMeasure : public DataTransferKit::ElementMeasure<
    DataTransferKit::MeshContainer<int> >
{
  public:

    MyMeasure( const DataTransferKit::MeshContainer<int>& mesh, 
		 const Teuchos::RCP< const Teuchos::Comm<int> >& comm )
	: d_mesh( mesh )
	, d_comm( comm )
    { /* ... */ }

    ~MyMeasure()
    { /* ... */ }

    // If the global id is valid, then set the element measure to 1, -1 if
    // invalid.
    Teuchos::Array<double> measure( 
	const Teuchos::ArrayRCP<
	    DataTransferKit::MeshContainer<int>::global_ordinal_type>& elements )
    {
	Teuchos::Array<double> measures( elements.size() );
	for ( int n = 0; n < elements.size(); ++n )
	{
	    if ( std::find( d_mesh.elementsBegin(),
			    d_mesh.elementsEnd(),
			    elements[n] ) != d_mesh.elementsEnd() )
	    {
		measures[n] = 1.0;
	    }
	    else
	    {
 		measures[n] = -1.0;
	    }
	}
	return measures;
    }

  private:

    DataTransferKit::MeshContainer<int>  d_mesh;
    Teuchos::RCP< const Teuchos::Comm<int> > d_comm;
};

//---------------------------------------------------------------------------//
// Mesh create functions.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> >  
buildHexMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some vertices.
    int num_vertices = edge_length*edge_length*2;
    int vertex_dim = 3;
    Teuchos::ArrayRCP<int> vertex_handles( num_vertices );
    Teuchos::ArrayRCP<double> coords( vertex_dim*num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    vertex_handles[ idx ] = (int) num_vertices*my_rank + idx + elem_offset;
	    coords[ idx ] = 2.0*my_rank + i/(edge_length-1);
	    coords[ num_vertices + idx ] = 2.0*my_rank + j/(edge_length-1);
	    coords[ 2*num_vertices + idx ] = 2.0*my_rank;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + num_vertices / 2;
	    vertex_handles[ idx ] = (int) num_vertices*my_rank + idx + elem_offset;
	    coords[ idx ] = 2.0*my_rank + i/(edge_length-1);
	    coords[ num_vertices + idx ] = 2.0*my_rank + j/(edge_length-1);
	    coords[ 2*num_vertices + idx ] = 2.0*my_rank + 1.0;
	}
    }

    // Make the hexahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::ArrayRCP<int> hex_handles( num_elements );
    Teuchos::ArrayRCP<int> hex_connectivity( 8*num_elements );
    int elem_idx, vertex_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    vertex_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    hex_handles[elem_idx] = num_elements*my_rank + elem_idx + elem_offset;

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
	new DataTransferKit::MeshContainer<int>( 3, vertex_handles, coords, 
						 DataTransferKit::DTK_HEXAHEDRON, 8,
						 hex_handles, hex_connectivity,
						 permutation_list ) );
}

//---------------------------------------------------------------------------//
// Geometry create functions. These geometries will span the entire domain,
// requiring them to be broadcast throughout the rendezvous.
//---------------------------------------------------------------------------//
void buildBoxGeometry( int my_size, int edge_size,
		       Teuchos::ArrayRCP<DataTransferKit::Box>& boxes,
		       Teuchos::ArrayRCP<int>& gids )
{
    Teuchos::ArrayRCP<DataTransferKit::Box> new_boxes(my_size);
    Teuchos::ArrayRCP<int> new_gids(my_size,0);
    for ( int i = 0; i < my_size; ++i )
    {
	new_gids[i] = i;
	double lo = 2.0*i;
	double hi =  2.0*i + 1.0;
	new_boxes[i] = DataTransferKit::Box( lo, lo, lo, hi, hi, hi );
    }
    boxes = new_boxes;
    gids = new_gids;
}

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntegralAssemblyMap, box_test )
{
    using namespace DataTransferKit;
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits<MeshType> MT;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Compute element ordinal offsets so we make unique global ordinals.
    int edge_size = 10;

    // Setup source mesh manager.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildHexMesh( my_rank, my_size, edge_size, 0 );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > source_mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 3 ) );

    // Setup target.
    int num_geom = my_size;
    int geometry_dim = 3;
    Teuchos::ArrayRCP<Box> geometry(0);
    Teuchos::ArrayRCP<int> geom_gids(0);
    int target_dim = 3;

    // Every proc will get the same geometry.
    buildBoxGeometry( my_size, edge_size, geometry, geom_gids );
    Teuchos::RCP<MyField> target_field = 
	Teuchos::rcp( new MyField( num_geom, target_dim ) );

    Teuchos::RCP< GeometryManager<Box,int> > target_geometry_manager =
	Teuchos::rcp( new GeometryManager<Box,int>( 
			  geometry, geom_gids, comm, geometry_dim ) );
    Teuchos::RCP<FieldManager<MyField> > target_space_manager = Teuchos::rcp( 
	new FieldManager<MyField>( target_field, comm ) );

    // Create field integrator and element measure.
    Teuchos::RCP< FieldIntegrator<MeshType ,MyField> > source_integrator;
    Teuchos::RCP<ElementMeasure<MeshType> > source_mesh_measure;
    source_integrator = Teuchos::rcp( new MyIntegrator( *mesh_blocks[0], comm ) );
    source_mesh_measure = Teuchos::rcp( new MyMeasure( *mesh_blocks[0], comm ) );

    // Create data target. This target is a scalar.
    // Setup and apply the integral assembly mapping.
    IntegralAssemblyMap<MeshType,Box> integral_assembly_map( 
	comm, source_mesh_manager->dim(), 1.0e-6, false );
    integral_assembly_map.setup( source_mesh_manager, source_mesh_measure,
				 target_geometry_manager );
    integral_assembly_map.apply( source_integrator, target_space_manager );

    // Check the integration. All elements in the mesh are in the box as this
    // is a true conformal situation and therefore the total integral should
    // be the global number of mesh elements. 
    for ( int d = 0; d < target_dim*num_geom; ++d )
    {
	TEST_EQUALITY( 2.0, target_field->getData()[d] );
    }
}

//---------------------------------------------------------------------------//
// end tstIntegralAssemblyMap6.cpp
//---------------------------------------------------------------------------//

