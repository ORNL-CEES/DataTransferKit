//---------------------------------------------------------------------------//
/*!
 * \file tstSharedDomainMap2.cpp
 * \author Stuart R. Slattery
 * \brief Shared domain map unit test 2.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_SharedDomainMap.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template <class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Mesh Implementation
//---------------------------------------------------------------------------//

class MyMesh
{
  public:
    typedef unsigned long int global_ordinal_type;

    MyMesh() { /* ... */}

    MyMesh( const Teuchos::Array<global_ordinal_type> &vertex_handles,
            const Teuchos::Array<double> &coords,
            const Teuchos::Array<global_ordinal_type> &element_handles,
            const Teuchos::Array<global_ordinal_type> &element_connectivity,
            const Teuchos::Array<int> &permutation_list )
        : d_vertex_handles( vertex_handles )
        , d_coords( coords )
        , d_element_handles( element_handles )
        , d_element_connectivity( element_connectivity )
        , d_permutation_list( permutation_list )
    { /* ... */
    }

    ~MyMesh() { /* ... */}

    Teuchos::Array<global_ordinal_type>::const_iterator verticesBegin() const
    {
        return d_vertex_handles.begin();
    }

    Teuchos::Array<global_ordinal_type>::const_iterator verticesEnd() const
    {
        return d_vertex_handles.end();
    }

    Teuchos::Array<double>::const_iterator coordsBegin() const
    {
        return d_coords.begin();
    }

    Teuchos::Array<double>::const_iterator coordsEnd() const
    {
        return d_coords.end();
    }

    Teuchos::Array<global_ordinal_type>::const_iterator elementsBegin() const
    {
        return d_element_handles.begin();
    }

    Teuchos::Array<global_ordinal_type>::const_iterator elementsEnd() const
    {
        return d_element_handles.end();
    }

    Teuchos::Array<global_ordinal_type>::const_iterator
    connectivityBegin() const
    {
        return d_element_connectivity.begin();
    }

    Teuchos::Array<global_ordinal_type>::const_iterator connectivityEnd() const
    {
        return d_element_connectivity.end();
    }

    Teuchos::Array<int>::const_iterator permutationBegin() const
    {
        return d_permutation_list.begin();
    }

    Teuchos::Array<int>::const_iterator permutationEnd() const
    {
        return d_permutation_list.end();
    }

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
        , d_data( dim * size )
    { /* ... */
    }

    ~MyField() { /* ... */}

    int dim() const { return d_dim; }

    size_type size() const { return d_data.size(); }

    bool empty() const { return d_data.empty(); }

    iterator begin() { return d_data.begin(); }

    const_iterator begin() const { return d_data.begin(); }

    iterator end() { return d_data.end(); }

    const_iterator end() const { return d_data.end(); }

    Teuchos::ArrayView<double> getData() { return d_data(); }

    const Teuchos::ArrayView<const double> getData() const { return d_data(); }

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
template <>
class MeshTraits<MyMesh>
{
  public:
    typedef MyMesh::global_ordinal_type global_ordinal_type;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator
        const_vertex_iterator;
    typedef Teuchos::Array<double>::const_iterator const_coordinate_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator
        const_element_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator
        const_connectivity_iterator;
    typedef Teuchos::Array<int>::const_iterator const_permutation_iterator;

    static inline int vertexDim( const MyMesh &mesh ) { return 3; }

    static inline const_vertex_iterator verticesBegin( const MyMesh &mesh )
    {
        return mesh.verticesBegin();
    }

    static inline const_vertex_iterator verticesEnd( const MyMesh &mesh )
    {
        return mesh.verticesEnd();
    }

    static inline const_coordinate_iterator coordsBegin( const MyMesh &mesh )
    {
        return mesh.coordsBegin();
    }

    static inline const_coordinate_iterator coordsEnd( const MyMesh &mesh )
    {
        return mesh.coordsEnd();
    }

    static inline DTK_ElementTopology elementTopology( const MyMesh &mesh )
    {
        return DTK_HEXAHEDRON;
    }

    static inline int verticesPerElement( const MyMesh &mesh ) { return 8; }

    static inline const_element_iterator elementsBegin( const MyMesh &mesh )
    {
        return mesh.elementsBegin();
    }

    static inline const_element_iterator elementsEnd( const MyMesh &mesh )
    {
        return mesh.elementsEnd();
    }

    static inline const_connectivity_iterator
    connectivityBegin( const MyMesh &mesh )
    {
        return mesh.connectivityBegin();
    }

    static inline const_connectivity_iterator
    connectivityEnd( const MyMesh &mesh )
    {
        return mesh.connectivityEnd();
    }

    static inline const_permutation_iterator
    permutationBegin( const MyMesh &mesh )
    {
        return mesh.permutationBegin();
    }

    static inline const_permutation_iterator
    permutationEnd( const MyMesh &mesh )
    {
        return mesh.permutationEnd();
    }
};

//---------------------------------------------------------------------------//
// Field Traits specification for MyField
template <>
class FieldTraits<MyField>
{
  public:
    typedef MyField field_type;
    typedef double value_type;
    typedef MyField::size_type size_type;
    typedef MyField::iterator iterator;
    typedef MyField::const_iterator const_iterator;

    static inline size_type dim( const MyField &field ) { return field.dim(); }

    static inline size_type size( const MyField &field )
    {
        return field.size();
    }

    static inline bool empty( const MyField &field ) { return field.empty(); }

    static inline iterator begin( MyField &field ) { return field.begin(); }

    static inline const_iterator begin( const MyField &field )
    {
        return field.begin();
    }

    static inline iterator end( MyField &field ) { return field.end(); }

    static inline const_iterator end( const MyField &field )
    {
        return field.end();
    }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// FieldEvaluator Implementation. This one populates a 3-vector.
class MyEvaluator
    : public DataTransferKit::FieldEvaluator<MyMesh::global_ordinal_type,
                                             MyField>
{
  public:
    MyEvaluator( const MyMesh &mesh,
                 const Teuchos::RCP<const Teuchos::Comm<int>> &comm )
        : d_mesh( mesh )
        , d_comm( comm )
    { /* ... */
    }

    ~MyEvaluator() { /* ... */}

    MyField
    evaluate( const Teuchos::ArrayRCP<MyMesh::global_ordinal_type> &elements,
              const Teuchos::ArrayRCP<double> &coords )
    {
        int num_elements = elements.size();
        MyField evaluated_data( num_elements, 3 );
        for ( int n = 0; n < num_elements; ++n )
        {
            if ( std::find( d_mesh.elementsBegin(), d_mesh.elementsEnd(),
                            elements[n] ) != d_mesh.elementsEnd() )
            {
                *( evaluated_data.begin() + n ) = d_comm->getRank() + 1.0;
                *( evaluated_data.begin() + num_elements + n ) =
                    d_comm->getRank() + 1.0;
                *( evaluated_data.begin() + 2 * num_elements + n ) =
                    d_comm->getRank() + 1.0;
            }
            else
            {
                *( evaluated_data.begin() + n ) = 0.0;
                *( evaluated_data.begin() + num_elements + n ) = 0.0;
                *( evaluated_data.begin() + 2 * num_elements + n ) = 0.0;
            }
        }
        return evaluated_data;
    }

  private:
    MyMesh d_mesh;
    Teuchos::RCP<const Teuchos::Comm<int>> d_comm;
};

//---------------------------------------------------------------------------//
// Mesh create functions.
//---------------------------------------------------------------------------//
Teuchos::RCP<MyMesh> buildMyMesh( int my_rank, int my_size, int edge_length )
{
    // Make some vertices.
    int num_vertices = edge_length * edge_length * 2;
    int vertex_dim = 3;
    Teuchos::Array<unsigned long int> vertex_handles( num_vertices );
    Teuchos::Array<double> coords( vertex_dim * num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
        for ( int i = 0; i < edge_length; ++i )
        {
            idx = i + j * edge_length;
            vertex_handles[idx] = (long int)num_vertices * my_rank + idx;
            coords[idx] = i + my_rank * ( edge_length - 1 );
            coords[num_vertices + idx] = j;
            coords[2 * num_vertices + idx] = 0.0;
        }
    }
    for ( int j = 0; j < edge_length; ++j )
    {
        for ( int i = 0; i < edge_length; ++i )
        {
            idx = i + j * edge_length + num_vertices / 2;
            vertex_handles[idx] = (long int)num_vertices * my_rank + idx;
            coords[idx] = i + my_rank * ( edge_length - 1 );
            coords[num_vertices + idx] = j;
            coords[2 * num_vertices + idx] = 1.0;
        }
    }

    // Make the hexahedrons.
    int num_elements = ( edge_length - 1 ) * ( edge_length - 1 );
    Teuchos::Array<unsigned long int> hex_handles( num_elements );
    Teuchos::Array<unsigned long int> hex_connectivity( 8 * num_elements );
    int elem_idx, vertex_idx;
    for ( int j = 0; j < ( edge_length - 1 ); ++j )
    {
        for ( int i = 0; i < ( edge_length - 1 ); ++i )
        {
            vertex_idx = i + j * edge_length;
            elem_idx = i + j * ( edge_length - 1 );

            hex_handles[elem_idx] = num_elements * my_rank + elem_idx;

            hex_connectivity[elem_idx] = vertex_handles[vertex_idx];

            hex_connectivity[num_elements + elem_idx] =
                vertex_handles[vertex_idx + 1];

            hex_connectivity[2 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + edge_length + 1];

            hex_connectivity[3 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + edge_length];

            hex_connectivity[4 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2];

            hex_connectivity[5 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2 + 1];

            hex_connectivity[6 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2 + edge_length + 1];

            hex_connectivity[7 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2 + edge_length];
        }
    }

    Teuchos::Array<int> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return Teuchos::rcp( new MyMesh( vertex_handles, coords, hex_handles,
                                     hex_connectivity, permutation_list ) );
}

//---------------------------------------------------------------------------//
// This function creates a tiled mesh such that the global mesh is not
// rectilinear.
Teuchos::RCP<MyMesh> buildTiledMesh( int my_rank, int my_size, int edge_length )
{
    // Make some vertices.
    int num_vertices = edge_length * edge_length * 2;
    int vertex_dim = 3;
    Teuchos::Array<unsigned long int> vertex_handles( num_vertices );
    Teuchos::Array<double> coords( vertex_dim * num_vertices );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
        for ( int i = 0; i < edge_length; ++i )
        {
            idx = i + j * edge_length;
            vertex_handles[idx] =
                (unsigned long int)num_vertices * my_rank + idx;
            coords[idx] = i + my_rank * ( edge_length - 1 );
            coords[num_vertices + idx] = j + my_rank * ( edge_length - 1 );
            coords[2 * num_vertices + idx] = 0.0;
        }
    }
    for ( int j = 0; j < edge_length; ++j )
    {
        for ( int i = 0; i < edge_length; ++i )
        {
            idx = i + j * edge_length + num_vertices / 2;
            vertex_handles[idx] =
                (unsigned long int)num_vertices * my_rank + idx;
            coords[idx] = i + my_rank * ( edge_length - 1 );
            coords[num_vertices + idx] = j + my_rank * ( edge_length - 1 );
            coords[2 * num_vertices + idx] = 1.0;
        }
    }

    // Make the hexahedrons.
    int num_elements = ( edge_length - 1 ) * ( edge_length - 1 );
    Teuchos::Array<unsigned long int> hex_handles( num_elements );
    Teuchos::Array<unsigned long int> hex_connectivity( 8 * num_elements );
    int elem_idx, vertex_idx;
    for ( int j = 0; j < ( edge_length - 1 ); ++j )
    {
        for ( int i = 0; i < ( edge_length - 1 ); ++i )
        {
            vertex_idx = i + j * edge_length;
            elem_idx = i + j * ( edge_length - 1 );

            hex_handles[elem_idx] = num_elements * my_rank + elem_idx;

            hex_connectivity[elem_idx] = vertex_handles[vertex_idx];

            hex_connectivity[num_elements + elem_idx] =
                vertex_handles[vertex_idx + 1];

            hex_connectivity[2 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + edge_length + 1];

            hex_connectivity[3 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + edge_length];

            hex_connectivity[4 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2];

            hex_connectivity[5 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2 + 1];

            hex_connectivity[6 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2 + edge_length + 1];

            hex_connectivity[7 * num_elements + elem_idx] =
                vertex_handles[vertex_idx + num_vertices / 2 + edge_length];
        }
    }

    Teuchos::Array<int> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
        permutation_list[i] = i;
    }

    return Teuchos::rcp( new MyMesh( vertex_handles, coords, hex_handles,
                                     hex_connectivity, permutation_list ) );
}

//---------------------------------------------------------------------------//
// Coordinate field create functions.
//---------------------------------------------------------------------------//
void buildCoordinateField( int my_rank, int my_size, int num_points,
                           int edge_size,
                           Teuchos::RCP<MyField> &coordinate_field )
{
    std::srand( my_rank * num_points * 2 );
    for ( int i = 0; i < num_points; ++i )
    {
        *( coordinate_field->begin() + i ) =
            my_size * ( edge_size - 1 ) * (double)std::rand() / RAND_MAX;
        *( coordinate_field->begin() + num_points + i ) =
            ( edge_size - 1 ) * (double)std::rand() / RAND_MAX;
        *( coordinate_field->begin() + 2 * num_points + i ) =
            (double)std::rand() / RAND_MAX;
    }
}

//---------------------------------------------------------------------------//
void buildExpandedCoordinateField( int my_rank, int my_size, int num_points,
                                   int edge_size,
                                   Teuchos::RCP<MyField> &coordinate_field )
{
    std::srand( my_rank * num_points * 2 );
    for ( int i = 0; i < num_points; ++i )
    {
        *( coordinate_field->begin() + i ) =
            my_size * ( edge_size ) * (double)std::rand() / RAND_MAX - 0.5;
        *( coordinate_field->begin() + num_points + i ) =
            ( edge_size ) * (double)std::rand() / RAND_MAX - 0.5;
        *( coordinate_field->begin() + 2 * num_points + i ) =
            1.2 * (double)std::rand() / RAND_MAX - 0.1;
    }
}

//---------------------------------------------------------------------------//
void buildTiledCoordinateField( int my_rank, int my_size, int num_points,
                                int edge_size,
                                Teuchos::RCP<MyField> &coordinate_field )
{
    std::srand( my_rank * num_points * 2 );
    for ( int i = 0; i < num_points; ++i )
    {
        *( coordinate_field->begin() + i ) =
            my_size * ( edge_size ) * (double)std::rand() / RAND_MAX - 0.5;
        *( coordinate_field->begin() + num_points + i ) =
            my_size * ( edge_size ) * (double)std::rand() / RAND_MAX - 0.5;
        *( coordinate_field->begin() + 2 * num_points + i ) =
            1.2 * (double)std::rand() / RAND_MAX - 0.1;
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
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MyMesh>> mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP<MeshManager<MyMesh>> source_mesh_manager =
        Teuchos::rcp( new MeshManager<MyMesh>( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field manager.
    int num_points = 1000;
    int point_dim = 3;
    Teuchos::RCP<MyField> coordinate_field =
        Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildCoordinateField( my_rank, my_size, num_points, edge_size,
                          coordinate_field );
    Teuchos::RCP<FieldManager<MyField>> target_coord_manager =
        Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<MyMesh::global_ordinal_type, MyField>>
        source_evaluator =
            Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );

    // Create data target. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field =
        Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP<FieldManager<MyField>> target_space_manager =
        Teuchos::rcp( new FieldManager<MyField>( target_field, comm ) );
    // Setup and apply the shared domain mapping.
    SharedDomainMap<MyMesh, MyField> shared_domain_map(
        comm, source_mesh_manager->dim() );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data.
    double source_rank;
    for ( int n = 0; n < num_points; ++n )
    {
        source_rank = std::floor( *( coordinate_field->begin() + n ) /
                                  ( edge_size - 1 ) );
        TEST_FLOATING_EQUALITY( source_rank + 1,
                                *( target_space_manager->field()->begin() + n ),
                                1.0e-14 );
        TEST_FLOATING_EQUALITY(
            source_rank + 1,
            *( target_space_manager->field()->begin() + n + num_points ),
            1.0e-14 );
        TEST_FLOATING_EQUALITY(
            source_rank + 1,
            *( target_space_manager->field()->begin() + n + 2 * num_points ),
            1.0e-14 );
    }
}

//---------------------------------------------------------------------------//
// Some points will be outside of the mesh in this test.
TEUCHOS_UNIT_TEST( SharedDomainMap, shared_domain_map_expanded_test2 )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MyMesh>> mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP<MeshManager<MyMesh>> source_mesh_manager =
        Teuchos::rcp( new MeshManager<MyMesh>( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field manager.
    int num_points = 1000;
    int point_dim = 3;
    Teuchos::RCP<MyField> coordinate_field =
        Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
                                  coordinate_field );
    Teuchos::RCP<FieldManager<MyField>> target_coord_manager =
        Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<MyMesh::global_ordinal_type, MyField>>
        source_evaluator =
            Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );

    // Create data target. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field =
        Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP<FieldManager<MyField>> target_space_manager =
        Teuchos::rcp( new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MyMesh, MyField> shared_domain_map(
        comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data if it is in the mesh and 0.0 if it is
    // outside.
    double source_rank;
    Teuchos::Array<unsigned long int> missing_points;
    for ( int n = 0; n < num_points; ++n )
    {
        if ( *( coordinate_field->begin() + n ) < 0.0 ||
             *( coordinate_field->begin() + n ) > ( edge_size - 1 ) * my_size ||
             *( coordinate_field->begin() + n + num_points ) < 0.0 ||
             *( coordinate_field->begin() + n + num_points ) > edge_size - 1 ||
             *( coordinate_field->begin() + n + 2 * num_points ) < 0.0 ||
             *( coordinate_field->begin() + n + 2 * num_points ) > 1.0 )
        {
            missing_points.push_back( n );
            TEST_EQUALITY( 0.0,
                           *( target_space_manager->field()->begin() + n ) );
            TEST_EQUALITY( 0.0, *( target_space_manager->field()->begin() + n +
                                   num_points ) );
            TEST_EQUALITY( 0.0, *( target_space_manager->field()->begin() + n +
                                   2 * num_points ) );
        }
        else
        {
            source_rank =
                std::floor( target_coord_manager->field()->getData()[n] /
                            ( edge_size - 1 ) );
            TEST_FLOATING_EQUALITY( source_rank + 1,
                                    target_space_manager->field()->getData()[n],
                                    1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager->field()->getData()[n + num_points],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager->field()->getData()[n + 2 * num_points],
                1.0e-14 );
        }
    }

    // Check the missing points.
    TEST_ASSERT( missing_points.size() > 0 );
    Teuchos::ArrayView<unsigned long int> missed_in_map =
        shared_domain_map.getMissedTargetPoints();
    TEST_EQUALITY( missing_points.size(), missed_in_map.size() );

    std::sort( missing_points.begin(), missing_points.end() );
    std::sort( missed_in_map.begin(), missed_in_map.end() );

    for ( int n = 0; n < (int)missing_points.size(); ++n )
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
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MyMesh>> mesh_blocks( 1 );
    mesh_blocks[0] = buildTiledMesh( my_rank, my_size, edge_size );
    Teuchos::RCP<MeshManager<MyMesh>> source_mesh_manager =
        Teuchos::rcp( new MeshManager<MyMesh>( mesh_blocks, comm, 3 ) );

    // Setup target coordinate field manager.
    int num_points = 1000;
    int point_dim = 3;
    Teuchos::RCP<MyField> coordinate_field =
        Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildTiledCoordinateField( my_rank, my_size, num_points, edge_size,
                               coordinate_field );
    Teuchos::RCP<FieldManager<MyField>> target_coord_manager =
        Teuchos::rcp( new FieldManager<MyField>( coordinate_field, comm ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<MyMesh::global_ordinal_type, MyField>>
        source_evaluator =
            Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );

    // Create data target. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field =
        Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP<FieldManager<MyField>> target_space_manager =
        Teuchos::rcp( new FieldManager<MyField>( target_field, comm ) );

    // Setup and apply the shared domain mapping.
    SharedDomainMap<MyMesh, MyField> shared_domain_map(
        comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager );
    shared_domain_map.apply( source_evaluator, target_space_manager );

    // Check the data transfer. Each target point should have been assigned
    // its source rank + 1 as data if it is in the mesh and 0.0 if it is
    // outside.
    double source_rank;
    Teuchos::Array<unsigned long int> missing_points;
    bool tagged;
    for ( int n = 0; n < num_points; ++n )
    {
        tagged = false;

        if ( *( coordinate_field->begin() + n ) < 0.0 ||
             *( coordinate_field->begin() + n ) > ( edge_size - 1 ) * my_size ||
             *( coordinate_field->begin() + n + num_points ) < 0.0 ||
             *( coordinate_field->begin() + n + num_points ) >
                 ( edge_size - 1 ) * my_size ||
             *( coordinate_field->begin() + n + 2 * num_points ) < 0.0 ||
             *( coordinate_field->begin() + n + 2 * num_points ) > 1.0 )
        {
            missing_points.push_back( n );
            TEST_EQUALITY( 0.0,
                           *( target_space_manager->field()->begin() + n ) );
            TEST_EQUALITY( 0.0, *( target_space_manager->field()->begin() + n +
                                   num_points ) );
            TEST_EQUALITY( 0.0, *( target_space_manager->field()->begin() + n +
                                   2 * num_points ) );
            tagged = true;
        }

        else
        {
            for ( int i = 0; i < my_size; ++i )
            {
                if ( *( coordinate_field->begin() + n ) >=
                         ( edge_size - 1 ) * i &&
                     *( coordinate_field->begin() + n ) <=
                         ( edge_size - 1 ) * ( i + 1 ) &&
                     *( coordinate_field->begin() + n + num_points ) >=
                         ( edge_size - 1 ) * i &&
                     *( coordinate_field->begin() + n + num_points ) <=
                         ( edge_size - 1 ) * ( i + 1 ) &&
                     *( coordinate_field->begin() + n + 2 * num_points ) >=
                         0.0 &&
                     *( coordinate_field->begin() + n + 2 * num_points ) <=
                         1.0 &&
                     !tagged )
                {
                    source_rank = std::floor(
                        target_coord_manager->field()->getData()[n] /
                        ( edge_size - 1 ) );
                    TEST_FLOATING_EQUALITY(
                        source_rank + 1,
                        target_space_manager->field()->getData()[n], 1.0e-14 );
                    TEST_FLOATING_EQUALITY( source_rank + 1,
                                            target_space_manager->field()
                                                ->getData()[n + num_points],
                                            1.0e-14 );
                    TEST_FLOATING_EQUALITY( source_rank + 1,
                                            target_space_manager->field()
                                                ->getData()[n + 2 * num_points],
                                            1.0e-14 );
                    tagged = true;
                }
            }

            if ( !tagged )
            {
                missing_points.push_back( n );
                TEST_EQUALITY(
                    0.0, *( target_space_manager->field()->begin() + n ) );
                TEST_EQUALITY( 0.0, *( target_space_manager->field()->begin() +
                                       n + num_points ) );
                TEST_EQUALITY( 0.0, *( target_space_manager->field()->begin() +
                                       n + 2 * num_points ) );
                tagged = true;
            }
        }

        TEST_ASSERT( tagged );
    }

    // Check the missing points.
    TEST_ASSERT( missing_points.size() > 0 );
    Teuchos::ArrayView<unsigned long int> missed_in_map =
        shared_domain_map.getMissedTargetPoints();
    TEST_EQUALITY( missing_points.size(), missed_in_map.size() );

    std::sort( missing_points.begin(), missing_points.end() );
    std::sort( missed_in_map.begin(), missed_in_map.end() );

    for ( int n = 0; n < (int)missing_points.size(); ++n )
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
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MyMesh>> mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP<MeshManager<MyMesh>> source_mesh_manager =
        Teuchos::rcp( new MeshManager<MyMesh>( mesh_blocks, comm, 3 ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<MyMesh::global_ordinal_type, MyField>>
        source_evaluator =
            Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );

    int num_points = 1000;
    int point_dim = 3;

    // Setup target coordinate field manager 1.
    Teuchos::RCP<MyField> coordinate_field_1 =
        Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
                                  coordinate_field_1 );
    Teuchos::RCP<FieldManager<MyField>> target_coord_manager_1 =
        Teuchos::rcp( new FieldManager<MyField>( coordinate_field_1, comm ) );

    // Create data target 1. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field_1 =
        Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP<FieldManager<MyField>> target_space_manager_1 =
        Teuchos::rcp( new FieldManager<MyField>( target_field_1, comm ) );

    // Setup target coordinate field manager 2.
    Teuchos::RCP<MyField> coordinate_field_2 =
        Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
                                  coordinate_field_2 );
    Teuchos::RCP<FieldManager<MyField>> target_coord_manager_2 =
        Teuchos::rcp( new FieldManager<MyField>( coordinate_field_2, comm ) );

    // Create data target 2. This target is a 3-vector.
    Teuchos::RCP<MyField> target_field_2 =
        Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP<FieldManager<MyField>> target_space_manager_2 =
        Teuchos::rcp( new FieldManager<MyField>( target_field_2, comm ) );

    // Setup and apply the shared domain mapping 1.
    SharedDomainMap<MyMesh, MyField> shared_domain_map_1(
        comm, source_mesh_manager->dim(), true );
    shared_domain_map_1.setup( source_mesh_manager, target_coord_manager_1 );
    shared_domain_map_1.apply( source_evaluator, target_space_manager_1 );

    // Check the first data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.
    double source_rank;
    Teuchos::Array<unsigned long int> missing_points_1;
    for ( int n = 0; n < num_points; ++n )
    {
        if ( *( coordinate_field_1->begin() + n ) < 0.0 ||
             *( coordinate_field_1->begin() + n ) >
                 ( edge_size - 1 ) * my_size ||
             *( coordinate_field_1->begin() + n + num_points ) < 0.0 ||
             *( coordinate_field_1->begin() + n + num_points ) >
                 edge_size - 1 ||
             *( coordinate_field_1->begin() + n + 2 * num_points ) < 0.0 ||
             *( coordinate_field_1->begin() + n + 2 * num_points ) > 1.0 )
        {
            missing_points_1.push_back( n );
            TEST_EQUALITY( 0.0,
                           *( target_space_manager_1->field()->begin() + n ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_1->field()->begin() +
                                   n + num_points ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_1->field()->begin() +
                                   n + 2 * num_points ) );
        }
        else
        {
            source_rank =
                std::floor( target_coord_manager_1->field()->getData()[n] /
                            ( edge_size - 1 ) );
            TEST_FLOATING_EQUALITY(
                source_rank + 1, target_space_manager_1->field()->getData()[n],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_1->field()->getData()[n + num_points],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_1->field()->getData()[n + 2 * num_points],
                1.0e-14 );
        }
    }

    // Check the missing points.
    TEST_ASSERT( missing_points_1.size() > 0 );
    Teuchos::ArrayView<unsigned long int> missed_in_map_1 =
        shared_domain_map_1.getMissedTargetPoints();
    TEST_EQUALITY( missing_points_1.size(), missed_in_map_1.size() );

    std::sort( missing_points_1.begin(), missing_points_1.end() );
    std::sort( missed_in_map_1.begin(), missed_in_map_1.end() );

    for ( int n = 0; n < (int)missing_points_1.size(); ++n )
    {
        TEST_EQUALITY( missing_points_1[n], missed_in_map_1[n] );
    }

    // Setup and apply the shared domain mapping 2.
    SharedDomainMap<MyMesh, MyField> shared_domain_map_2(
        comm, source_mesh_manager->dim(), true );
    shared_domain_map_2.setup( source_mesh_manager, target_coord_manager_2 );
    shared_domain_map_2.apply( source_evaluator, target_space_manager_2 );

    // Check the second data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.
    Teuchos::Array<unsigned long int> missing_points_2;
    for ( int n = 0; n < num_points; ++n )
    {
        if ( *( coordinate_field_2->begin() + n ) < 0.0 ||
             *( coordinate_field_2->begin() + n ) >
                 ( edge_size - 1 ) * my_size ||
             *( coordinate_field_2->begin() + n + num_points ) < 0.0 ||
             *( coordinate_field_2->begin() + n + num_points ) >
                 edge_size - 1 ||
             *( coordinate_field_2->begin() + n + 2 * num_points ) < 0.0 ||
             *( coordinate_field_2->begin() + n + 2 * num_points ) > 1.0 )
        {
            missing_points_2.push_back( n );
            TEST_EQUALITY( 0.0,
                           *( target_space_manager_2->field()->begin() + n ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_2->field()->begin() +
                                   n + num_points ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_2->field()->begin() +
                                   n + 2 * num_points ) );
        }
        else
        {
            source_rank =
                std::floor( target_coord_manager_2->field()->getData()[n] /
                            ( edge_size - 1 ) );
            TEST_FLOATING_EQUALITY(
                source_rank + 1, target_space_manager_2->field()->getData()[n],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_2->field()->getData()[n + num_points],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_2->field()->getData()[n + 2 * num_points],
                1.0e-14 );
        }
    }

    // Check the missing points 2.
    TEST_ASSERT( missing_points_2.size() > 0 );
    Teuchos::ArrayView<unsigned long int> missed_in_map_2 =
        shared_domain_map_2.getMissedTargetPoints();
    TEST_EQUALITY( missing_points_2.size(), missed_in_map_2.size() );

    std::sort( missing_points_2.begin(), missing_points_2.end() );
    std::sort( missed_in_map_2.begin(), missed_in_map_2.end() );

    for ( int n = 0; n < (int)missing_points_2.size(); ++n )
    {
        TEST_EQUALITY( missing_points_2[n], missed_in_map_2[n] );
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
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh manager.
    int edge_size = 4;
    Teuchos::ArrayRCP<Teuchos::RCP<MyMesh>> mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh( my_rank, my_size, edge_size );
    Teuchos::RCP<MeshManager<MyMesh>> source_mesh_manager =
        Teuchos::rcp( new MeshManager<MyMesh>( mesh_blocks, comm, 3 ) );

    // Create field evaluator.
    Teuchos::RCP<FieldEvaluator<MyMesh::global_ordinal_type, MyField>>
        source_evaluator =
            Teuchos::rcp( new MyEvaluator( *mesh_blocks[0], comm ) );

    int num_points = 1000;
    int point_dim = 3;

    // Setup target coordinate field manager 1.
    Teuchos::RCP<MyField> coordinate_field_1 =
        Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
                                  coordinate_field_1 );
    Teuchos::RCP<FieldManager<MyField>> target_coord_manager_1 =
        Teuchos::rcp( new FieldManager<MyField>( coordinate_field_1, comm ) );

    // Create data target 1. This target is a 3-vector.
    int target_dim = 3;
    Teuchos::RCP<MyField> target_field_1 =
        Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP<FieldManager<MyField>> target_space_manager_1 =
        Teuchos::rcp( new FieldManager<MyField>( target_field_1, comm ) );

    // Setup target coordinate field manager 2.
    Teuchos::RCP<MyField> coordinate_field_2 =
        Teuchos::rcp( new MyField( num_points, point_dim ) );
    buildExpandedCoordinateField( my_rank, my_size, num_points, edge_size,
                                  coordinate_field_2 );
    Teuchos::RCP<FieldManager<MyField>> target_coord_manager_2 =
        Teuchos::rcp( new FieldManager<MyField>( coordinate_field_2, comm ) );

    // Create data target 2. This target is a 3-vector.
    Teuchos::RCP<MyField> target_field_2 =
        Teuchos::rcp( new MyField( num_points, target_dim ) );
    Teuchos::RCP<FieldManager<MyField>> target_space_manager_2 =
        Teuchos::rcp( new FieldManager<MyField>( target_field_2, comm ) );

    // Setup shared domain mapping.
    SharedDomainMap<MyMesh, MyField> shared_domain_map(
        comm, source_mesh_manager->dim(), true );
    shared_domain_map.setup( source_mesh_manager, target_coord_manager_1 );

    // Apply the shared domain mapping 1.
    shared_domain_map.apply( source_evaluator, target_space_manager_1 );

    // Check the first data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.
    double source_rank;
    Teuchos::Array<unsigned long int> missing_points_1;
    for ( int n = 0; n < num_points; ++n )
    {
        if ( *( coordinate_field_1->begin() + n ) < 0.0 ||
             *( coordinate_field_1->begin() + n ) >
                 ( edge_size - 1 ) * my_size ||
             *( coordinate_field_1->begin() + n + num_points ) < 0.0 ||
             *( coordinate_field_1->begin() + n + num_points ) >
                 edge_size - 1 ||
             *( coordinate_field_1->begin() + n + 2 * num_points ) < 0.0 ||
             *( coordinate_field_1->begin() + n + 2 * num_points ) > 1.0 )
        {
            missing_points_1.push_back( n );
            TEST_EQUALITY( 0.0,
                           *( target_space_manager_1->field()->begin() + n ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_1->field()->begin() +
                                   n + num_points ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_1->field()->begin() +
                                   n + 2 * num_points ) );
        }
        else
        {
            source_rank =
                std::floor( target_coord_manager_1->field()->getData()[n] /
                            ( edge_size - 1 ) );
            TEST_FLOATING_EQUALITY(
                source_rank + 1, target_space_manager_1->field()->getData()[n],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_1->field()->getData()[n + num_points],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_1->field()->getData()[n + 2 * num_points],
                1.0e-14 );
        }
    }

    // Check the missing points.
    TEST_ASSERT( missing_points_1.size() > 0 );
    Teuchos::ArrayView<unsigned long int> missed_in_map_1 =
        shared_domain_map.getMissedTargetPoints();
    TEST_EQUALITY( missing_points_1.size(), missed_in_map_1.size() );

    std::sort( missing_points_1.begin(), missing_points_1.end() );
    std::sort( missed_in_map_1.begin(), missed_in_map_1.end() );

    for ( int n = 0; n < (int)missing_points_1.size(); ++n )
    {
        TEST_EQUALITY( missing_points_1[n], missed_in_map_1[n] );
    }

    // Apply the shared domain mapping 2.
    shared_domain_map.apply( source_evaluator, target_space_manager_2 );

    // Check the second data transfer. Each target point should have been
    // assigned its source rank + 1 as data if it is in the mesh and 0.0 if it
    // is outside.
    Teuchos::Array<unsigned long int> missing_points_2;
    for ( int n = 0; n < num_points; ++n )
    {
        if ( *( coordinate_field_2->begin() + n ) < 0.0 ||
             *( coordinate_field_2->begin() + n ) >
                 ( edge_size - 1 ) * my_size ||
             *( coordinate_field_2->begin() + n + num_points ) < 0.0 ||
             *( coordinate_field_2->begin() + n + num_points ) >
                 edge_size - 1 ||
             *( coordinate_field_2->begin() + n + 2 * num_points ) < 0.0 ||
             *( coordinate_field_2->begin() + n + 2 * num_points ) > 1.0 )
        {
            missing_points_2.push_back( n );
            TEST_EQUALITY( 0.0,
                           *( target_space_manager_2->field()->begin() + n ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_2->field()->begin() +
                                   n + num_points ) );
            TEST_EQUALITY( 0.0, *( target_space_manager_2->field()->begin() +
                                   n + 2 * num_points ) );
        }
        else
        {
            source_rank =
                std::floor( target_coord_manager_2->field()->getData()[n] /
                            ( edge_size - 1 ) );
            TEST_FLOATING_EQUALITY(
                source_rank + 1, target_space_manager_2->field()->getData()[n],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_2->field()->getData()[n + num_points],
                1.0e-14 );
            TEST_FLOATING_EQUALITY(
                source_rank + 1,
                target_space_manager_2->field()->getData()[n + 2 * num_points],
                1.0e-14 );
        }
    }

    // Check the missing points 2.
    TEST_ASSERT( missing_points_2.size() > 0 );
    Teuchos::ArrayView<unsigned long int> missed_in_map_2 =
        shared_domain_map.getMissedTargetPoints();
    TEST_EQUALITY( missing_points_2.size(), missed_in_map_2.size() );

    std::sort( missing_points_2.begin(), missing_points_2.end() );
    std::sort( missed_in_map_2.begin(), missed_in_map_2.end() );

    for ( int n = 0; n < (int)missing_points_2.size(); ++n )
    {
        TEST_EQUALITY( missing_points_2[n], missed_in_map_2[n] );
    }
}

//---------------------------------------------------------------------------//
// end tstSharedDomainMap2.cpp
//---------------------------------------------------------------------------//
