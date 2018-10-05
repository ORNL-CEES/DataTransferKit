// The C API tests try to follow tests in tstUserApplication
#include "DTK_APIConstants.h"
#include <DTK_C_API.h>

#include <mpi.h>

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//---------------------------------------------------------------------------//
// User class
//---------------------------------------------------------------------------//
typedef struct UserTestClass
{
    unsigned _space_dim;
    unsigned _size_1;
    unsigned _size_2;
    unsigned _offset;
    const char *_field_name;
    double *_data;
} UserTestClass;

//---------------------------------------------------------------------------//
// User functions.
//---------------------------------------------------------------------------//
// Get the size parameters for building a node list.
void node_list_size( void *user_data, unsigned *space_dim,
                     size_t *local_num_nodes )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *space_dim = u->_space_dim;
    *local_num_nodes = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a node list.
void node_list_data( void *user_data, Coordinate *coordinates )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( int i = 0; i < u->_size_1; ++i )
    {
        for ( int j = 0; j < u->_space_dim; j++ )
            coordinates[j * u->_size_1 + i] = i + j + u->_offset;
    }
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a bounding volume list.
void bounding_volume_list_size( void *user_data, unsigned *space_dim,
                                size_t *local_num_volumes )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *space_dim = u->_space_dim;
    *local_num_volumes = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a bounding volume list.
void bounding_volume_list_data( void *user_data, Coordinate *bounding_volumes )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( size_t v = 0; v < u->_size_1; v++ )
    {
        for ( unsigned d = 0; d < u->_space_dim; ++d )
        {
            for ( unsigned h = 0; h < 2; ++h )
            {
                unsigned index =
                    u->_size_1 * u->_space_dim * h + u->_size_1 * d + v;
                bounding_volumes[index] = v + d + h + u->_offset;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a polyhedron list.
void polyhedron_list_size( void *user_data, unsigned *space_dim,
                           size_t *local_num_nodes, size_t *local_num_faces,
                           size_t *total_face_nodes, size_t *local_num_cells,
                           size_t *total_cell_faces )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *space_dim = u->_space_dim;
    *local_num_nodes = u->_size_1;
    *local_num_faces = u->_size_1;
    *total_face_nodes = u->_size_1;
    *local_num_cells = u->_size_1;
    *total_cell_faces = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a polyhedron list.
void polyhedron_list_data( void *user_data, Coordinate *coordinates,
                           LocalOrdinal *faces, unsigned *nodes_per_face,
                           LocalOrdinal *cells, unsigned *faces_per_cell,
                           int *face_orientation )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        for ( unsigned d = 0; d < u->_space_dim; ++d )
        {
            coordinates[u->_size_1 * d + n] = n + d + u->_offset;
        }
        faces[n] = n + u->_offset;
        nodes_per_face[n] = n + u->_offset;
        cells[n] = n + u->_offset;
        faces_per_cell[n] = n + u->_offset;
        face_orientation[n] = 1;
    }
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a cell list.
void cell_list_size( void *user_data, unsigned *space_dim,
                     size_t *local_num_nodes, size_t *local_num_cells,
                     size_t *total_cell_nodes )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *space_dim = u->_space_dim;
    *local_num_nodes = u->_size_1;
    *local_num_cells = u->_size_1;
    *total_cell_nodes = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a cell list.
void cell_list_data( void *user_data, Coordinate *coordinates,
                     LocalOrdinal *cells, DTK_CellTopology *cell_topologies )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        for ( unsigned d = 0; d < u->_space_dim; ++d )
            coordinates[u->_size_1 * d + n] = n + d + u->_offset;
        cells[n] = n + u->_offset;
        cell_topologies[n] = DTK_TET_4;
    }
}

//---------------------------------------------------------------------------//
// Get the size parameters for a boundary.
void boundary_size( void *user_data, size_t *local_num_faces )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *local_num_faces = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a boundary.
void boundary_data( void *user_data, LocalOrdinal *boundary_cells,
                    unsigned *cell_faces_on_boundary )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        boundary_cells[n] = n + u->_offset;
        cell_faces_on_boundary[n] = n + u->_offset;
    }
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a cell list.
void adjacency_list_size( void *user_data, size_t *total_adjacencies )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *total_adjacencies = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a adjacency list.
void adjacency_list_data( void *user_data, GlobalOrdinal *global_cell_ids,
                          GlobalOrdinal *adjacent_global_cell_ids,
                          unsigned *adjacencies_per_cell )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        global_cell_ids[n] = n + u->_offset;
        adjacent_global_cell_ids[n] = n;
        adjacencies_per_cell[n] = 1;
    }
}

//---------------------------------------------------------------------------//
// Get the size parameters for a degree-of-freedom id map with a single
// number of dofs per object.
void dof_map_size( void *user_data, size_t *local_num_dofs,
                   size_t *local_num_objects, unsigned *dofs_per_object )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *local_num_dofs = u->_size_1;
    *local_num_objects = u->_size_1;
    *dofs_per_object = u->_size_2;
}

//---------------------------------------------------------------------------//
// Get the data for a degree-of-freedom id map with a single number of
// dofs per object.
void dof_map_data( void *user_data, GlobalOrdinal *global_dof_ids,
                   LocalOrdinal *object_dof_ids, char *discretization_type )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        global_dof_ids[n] = n + u->_offset;
        for ( unsigned d = 0; d < u->_size_2; ++d )
            object_dof_ids[u->_size_1 * d + n] = n + d + u->_offset;
    }

    strcpy( discretization_type, "unit_test_discretization" );
}

//---------------------------------------------------------------------------//
// Get the size parameters for a degree-of-freedom id map with a
// multiple number of dofs per object (e.g. mixed topology cell lists or
// polyhedron lists).
void mixed_topology_dof_map_size( void *user_data, size_t *local_num_dofs,
                                  size_t *local_num_objects,
                                  size_t *total_dofs_per_object )
{
    UserTestClass *u = (UserTestClass *)user_data;

    *local_num_dofs = u->_size_1;
    *local_num_objects = u->_size_1;
    *total_dofs_per_object = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a multiple object degree-of-freedom id map
// (e.g. mixed topology cell lists or polyhedron lists).
void mixed_topology_dof_map_data( void *user_data,
                                  GlobalOrdinal *global_dof_ids,
                                  LocalOrdinal *object_dof_ids,
                                  unsigned *dofs_per_object,
                                  char *discretization_type )
{
    UserTestClass *u = (UserTestClass *)user_data;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        global_dof_ids[n] = n + u->_offset;
        object_dof_ids[n] = n + u->_offset;
        dofs_per_object[n] = u->_size_2;
    }

    strcpy( discretization_type, "unit_test_discretization" );
}

//---------------------------------------------------------------------------//
// Get the size parameters for a field. Field must be of size
// local_num_dofs in the associated dof_id_map.
void field_size( void *user_data, const char *field_name,
                 unsigned *field_dimension, size_t *local_num_dofs )
{
    UserTestClass *u = (UserTestClass *)user_data;

    // Here one could do actions depending on the name, but in the tests we
    // simply ignore it
    (void)field_name;

    *field_dimension = u->_space_dim;
    *local_num_dofs = u->_size_1;
}

//---------------------------------------------------------------------------//
// Pull data from application into a field.
void pull_field_data( void *user_data, const char *field_name,
                      double *field_dofs )

{
    UserTestClass *u = (UserTestClass *)user_data;

    // Here one could do actions depending on the name, but in the tests we
    // simply ignore it
    (void)field_name;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        for ( unsigned d = 0; d < u->_space_dim; ++d )
            field_dofs[d * u->_size_1 + n] = u->_data[d * u->_size_1 + n];
    }
}

//---------------------------------------------------------------------------//
// Push data from a field into the application.
void push_field_data( void *user_data, const char *field_name,
                      const double *field_dofs )

{
    UserTestClass *u = (UserTestClass *)user_data;

    // Here one could do actions depending on the name, but in the tests we
    // simply ignore it
    (void)field_name;

    for ( size_t n = 0; n < u->_size_1; n++ )
    {
        for ( unsigned d = 0; d < u->_space_dim; ++d )
            u->_data[d * u->_size_1 + n] = field_dofs[d * u->_size_1 + n];
    }
}

//---------------------------------------------------------------------------//
// Evaluate a field at a given set of points in a given set of objects.
void evaluate_field( void *user_data, const char *field_name,
                     const size_t num_points, const unsigned space_dim,
                     const Coordinate *evaluation_points,
                     const LocalOrdinal *object_ids, double *values )
{
    // Here one could do actions depending on the name, but in the tests we
    // simply ignore it
    (void)field_name;

    for ( size_t n = 0; n < num_points; n++ )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            values[d * num_points + n] =
                evaluation_points[d * num_points + n] + object_ids[n];
    }
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
extern int check_registry( const char *test_name,
                           DTK_UserApplicationHandle handle );

int test_node_list( DTK_UserApplicationHandle dtk_handle, UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_NODE_LIST_SIZE_FUNCTION,
                         node_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_NODE_LIST_DATA_FUNCTION,
                         node_list_data, &u );

    return check_registry( "test_node_list", dtk_handle );
}

int test_bounding_volume_list( DTK_UserApplicationHandle dtk_handle,
                               UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_BOUNDING_VOLUME_LIST_SIZE_FUNCTION,
                         bounding_volume_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_BOUNDING_VOLUME_LIST_DATA_FUNCTION,
                         bounding_volume_list_data, &u );

    return check_registry( "test_bounding_volume_list", dtk_handle );
}

int test_polyhedron_list( DTK_UserApplicationHandle dtk_handle,
                          UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_POLYHEDRON_LIST_SIZE_FUNCTION,
                         polyhedron_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_POLYHEDRON_LIST_DATA_FUNCTION,
                         polyhedron_list_data, &u );

    return check_registry( "test_polyhedron_list", dtk_handle );
}

int test_multiple_topology_cell( DTK_UserApplicationHandle dtk_handle,
                                 UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_CELL_LIST_SIZE_FUNCTION,
                         cell_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_CELL_LIST_DATA_FUNCTION,
                         cell_list_data, &u );

    return check_registry( "test_multiple_topology_cell", dtk_handle );
}

int test_boundary( DTK_UserApplicationHandle dtk_handle, UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_BOUNDARY_SIZE_FUNCTION, boundary_size,
                         &u );
    DTK_setUserFunction( dtk_handle, DTK_BOUNDARY_DATA_FUNCTION, boundary_data,
                         &u );
    DTK_setUserFunction( dtk_handle, DTK_CELL_LIST_SIZE_FUNCTION,
                         cell_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_CELL_LIST_DATA_FUNCTION,
                         cell_list_data, &u );
    DTK_setUserFunction( dtk_handle, DTK_POLYHEDRON_LIST_SIZE_FUNCTION,
                         polyhedron_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_POLYHEDRON_LIST_DATA_FUNCTION,
                         polyhedron_list_data, &u );

    return check_registry( "test_boundary", dtk_handle );
}

int test_adjacency_list( DTK_UserApplicationHandle dtk_handle, UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_ADJACENCY_LIST_SIZE_FUNCTION,
                         adjacency_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_ADJACENCY_LIST_DATA_FUNCTION,
                         adjacency_list_data, &u );
    DTK_setUserFunction( dtk_handle, DTK_CELL_LIST_SIZE_FUNCTION,
                         cell_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_CELL_LIST_DATA_FUNCTION,
                         cell_list_data, &u );
    DTK_setUserFunction( dtk_handle, DTK_POLYHEDRON_LIST_SIZE_FUNCTION,
                         polyhedron_list_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_POLYHEDRON_LIST_DATA_FUNCTION,
                         polyhedron_list_data, &u );

    return check_registry( "test_adjacency_list", dtk_handle );
}

int test_single_topology_dof( DTK_UserApplicationHandle dtk_handle,
                              UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_DOF_MAP_SIZE_FUNCTION, dof_map_size,
                         &u );
    DTK_setUserFunction( dtk_handle, DTK_DOF_MAP_DATA_FUNCTION, dof_map_data,
                         &u );

    return check_registry( "test_single_topology_dof", dtk_handle );
}

int test_multiple_topology_dof( DTK_UserApplicationHandle dtk_handle,
                                UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION,
                         mixed_topology_dof_map_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION,
                         mixed_topology_dof_map_data, &u );

    return check_registry( "test_multiple_topology_dof", dtk_handle );
}

int test_field_push_pull( DTK_UserApplicationHandle dtk_handle,
                          UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_FIELD_SIZE_FUNCTION, field_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_PULL_FIELD_DATA_FUNCTION,
                         pull_field_data, &u );
    DTK_setUserFunction( dtk_handle, DTK_PUSH_FIELD_DATA_FUNCTION,
                         push_field_data, &u );

    return check_registry( "test_field_push_pull", dtk_handle );
}

int test_field_eval( DTK_UserApplicationHandle dtk_handle, UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_FIELD_SIZE_FUNCTION, field_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_EVALUATE_FIELD_FUNCTION,
                         evaluate_field, &u );

    return check_registry( "test_field_eval", dtk_handle );
}

int test_missing_function( DTK_UserApplicationHandle dtk_handle,
                           UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_NODE_LIST_SIZE_FUNCTION,
                         node_list_size, &u );

    return check_registry( "test_missing_function", dtk_handle );
}

int test_too_many_functions( DTK_UserApplicationHandle dtk_handle,
                             UserTestClass u )
{
    DTK_setUserFunction( dtk_handle, DTK_DOF_MAP_SIZE_FUNCTION, dof_map_size,
                         &u );
    DTK_setUserFunction( dtk_handle, DTK_DOF_MAP_DATA_FUNCTION, dof_map_data,
                         &u );
    DTK_setUserFunction( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION,
                         mixed_topology_dof_map_size, &u );
    DTK_setUserFunction( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION,
                         mixed_topology_dof_map_data, &u );

    return check_registry( "test_too_many_functions", dtk_handle );
}

int main( int argc, char *argv[] )
{
    MPI_Init( &argc, &argv );

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    DTK_MemorySpace memory_space = DTK_HOST_SPACE;

    int opt;
    while ( ( opt = getopt( argc, argv, "s:h" ) ) != -1 )
    {
        switch ( opt )
        {
        case 's':
            if ( !strcmp( optarg, "serial" ) )
                memory_space = DTK_HOST_SPACE;
            else if ( !strcmp( optarg, "openmp" ) )
                memory_space = DTK_HOST_SPACE;
            else if ( !strcmp( optarg, "cuda" ) )
                memory_space = DTK_CUDAUVM_SPACE;
            else
            {
                printf( "Unknown memory space\n" );
                return EXIT_FAILURE;
            }
            break;

        case 'h':
            printf( "Usage: %s [-s <serial|openmp|cuda>] [-h]\n", argv[0] );
            return EXIT_FAILURE;
        }
    }

    if ( !comm_rank )
    {
        printf( "DTK version: %s\n", DTK_version() );
        printf( "DTK hash: %s\n", DTK_gitCommitHash() );
    }

    UserTestClass u;
    u._space_dim = SPACE_DIM;
    u._size_1 = SIZE_1;
    u._size_2 = SIZE_2;
    u._offset = OFFSET;
    u._field_name = FIELD_NAME;
    u._data = (double *)calloc( u._size_1 * u._space_dim, sizeof( double ) );

    int rv = 0;

    {
        DTK_createUserApplication( memory_space );
        rv |= ( errno == 0 );
        const char *errormsg = DTK_error( errno );
        rv |= ( strcmp( errormsg, "DTK error: DTK is not initialized" ) );
    }

    DTK_initializeCmd( &argc, &argv );
    rv |= ( DTK_isInitialized() ? 0 : 1 );

    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= ( errno != 0 );
        const char *errormsg = DTK_error( errno );
        rv |= ( strcmp( errormsg, "" ) );

        rv |= ( DTK_isValidUserApplication( dtk_handle ) ? 0 : 1 );
        DTK_destroyUserApplication( dtk_handle );
        rv |= ( DTK_isValidUserApplication( dtk_handle ) ? 1 : 0 );
    }
    {
        DTK_UserApplicationHandle dtk_handle;
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wuninitialized"
#else
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
        rv |= ( DTK_isValidUserApplication( dtk_handle ) ? 1 : 0 );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle;
        DTK_setUserFunction( dtk_handle, DTK_NODE_LIST_SIZE_FUNCTION,
                             node_list_size, &u );
        rv |= ( errno == 0 );

        const char *errormsg = DTK_error( errno );
        rv |= ( strcmp( errormsg, "DTK error: invalid DTK handle" ) );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_node_list( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_bounding_volume_list( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_polyhedron_list( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_multiple_topology_cell( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_boundary( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_adjacency_list( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_single_topology_dof( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_multiple_topology_dof( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_field_push_pull( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_field_eval( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_missing_function( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }
    {
        DTK_UserApplicationHandle dtk_handle =
            DTK_createUserApplication( memory_space );
        rv |= test_too_many_functions( dtk_handle, u );
        DTK_destroyUserApplication( dtk_handle );
    }

    DTK_finalize();

    if ( !comm_rank )
    {
        if ( rv == 0 )
            printf( "End Result: TEST PASSED\n" );
        else
            printf( "End Result: TEST FAILED\n" );
    }

    free( u._data );

    MPI_Finalize();

    return EXIT_SUCCESS;
}
