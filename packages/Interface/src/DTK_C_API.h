/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file
 * \brief C interface header.
 */
#ifndef DTK_C_API_H
#define DTK_C_API_H

#include <DataTransferKit_config.hpp>

#include <DTK_Types.h>

#ifndef __cplusplus
#include <stdbool.h> // for bool
#include <stddef.h>  // for size_t
#else
#include <cstddef> // for size_t
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! Return current DTK version. */
extern const char *DTK_version();

/*! Get current repository hash.
 *
 * Returns an error string if it is not a GIT repository
 */
extern const char *DTK_git_commit_hash();

// Use incomplete types to differentiate between handles.
// We never define the incomplete structs.
typedef struct _DTK_UserApplicationHandle *DTK_UserApplicationHandle;

typedef enum { DTK_SERIAL, DTK_OPENMP, DTK_CUDA } DTK_ExecutionSpace;

/*! Create a DTK handle. */
extern DTK_UserApplicationHandle DTK_create( DTK_ExecutionSpace space );

/*! Check if the DTK handle is valid */
extern bool DTK_is_valid( DTK_UserApplicationHandle handle );

/*! Destroy a DTK handle */
extern void DTK_destroy( DTK_UserApplicationHandle handle );

/**
 * \defgroup c_interface_to_dtk_core Initialize/finalize DTK
 * @{
 */

/** \brief Initializes the DTK execution environment.
 *
 *  This initializes Kokkos if it has not already been initialized.
 */
extern void DTK_initialize();
/** \brief Initialize DTK.
 *
 *  This initializes Kokkos if it has not already been initialized.
 *
 *  \param argc Pointer to the number of argument.
 *  \param argv Pointer to the argument vector.
 *
 *  This version of initialize() effectively calls <code>Kokkos::initialize(
 *  *argc, *argv )</code>.  Pointers to \c argc and \c argv arguments are passed
 *  in order to match MPI_Init's interface.  This function name was suffixed with
 *  _cmd because, unlike C++, C does not allow to overload functions.
 *
 *  <!--
 * FIXME: provide a version based on full Kokkos::InitArguments
 * FIXME: can we do something with the name here?
 *  -->
 */
extern void DTK_initialize_cmd( int *argc, char ***argv );

/** \brief Indicates whether DTK has been initialized.
 *
 *  This function may be used to determine whether DTK has been initialized.
 */
extern bool DTK_is_initialized();

/** \brief Finalizes DTK.
 *
 *  This function terminates the DTK execution environment.  If DTK initialized
 *  Kokkos, finalize Kokkos.  However, if Kokkos was initialized before DTK,
 *  then this function does NOT finalize Kokkos.
 */
extern void DTK_finalize();

/**@}*/

typedef enum {
    DTK_NODE_LIST_SIZE_FUNCTION,
    DTK_NODE_LIST_DATA_FUNCTION,
    DTK_BOUNDING_VOLUME_LIST_SIZE_FUNCTION,
    DTK_BOUNDING_VOLUME_LIST_DATA_FUNCTION,
    DTK_POLYHEDRON_LIST_SIZE_FUNCTION,
    DTK_POLYHEDRON_LIST_DATA_FUNCTION,
    DTK_CELL_LIST_SIZE_FUNCTION,
    DTK_CELL_LIST_DATA_FUNCTION,
    DTK_MIXED_TOPOLOGY_CELL_LIST_SIZE_FUNCTION,
    DTK_MIXED_TOPOLOGY_CELL_LIST_DATA_FUNCTION,
    DTK_BOUNDARY_SIZE_FUNCTION,
    DTK_BOUNDARY_DATA_FUNCTION,
    DTK_DOF_MAP_SIZE_FUNCTION,
    DTK_DOF_MAP_DATA_FUNCTION,
    DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION,
    DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION,
    DTK_FIELD_SIZE_FUNCTION,
    DTK_PULL_FIELD_DATA_FUNCTION,
    DTK_PUSH_FIELD_DATA_FUNCTION,
    DTK_EVALUATE_FIELD_FUNCTION
} DTK_FunctionType;

/*! Register a function with DTK. */
extern void DTK_set_function( DTK_UserApplicationHandle handle,
                              DTK_FunctionType type, void ( *f )(),
                              void *user_data );

/*! Get the size parameters for building a node list. */
typedef void ( *DTK_NodeListSizeFunction )( void *user_data,
                                            unsigned *space_dim,
                                            size_t *local_num_nodes,
                                            bool *has_ghosts );

/*! Get the data for a node list. */
typedef void ( *DTK_NodeListDataFunction )( void *user_data,
                                            Coordinate *coordinates,
                                            bool *is_ghost_node );

/*! Get the size parameters for building a bounding volume list. */
typedef void ( *DTK_BoundingVolumeListSizeFunction )( void *user_data,
                                                      unsigned *space_dim,
                                                      size_t *local_num_volumes,
                                                      bool *has_ghosts );

/*! Get the data for a bounding volume list. */
typedef void ( *DTK_BoundingVolumeListDataFunction )(
    void *user_data, Coordinate *bounding_volumes, bool *is_ghost_volume );

/*! Get the size parameters for building a polyhedron list. */
typedef void ( *DTK_PolyhedronListSizeFunction )(
    void *user_data, unsigned *space_dim, size_t *local_num_nodes,
    size_t *local_num_faces, size_t *total_nodes_per_face,
    size_t *local_num_cells, size_t *total_faces_per_cell, bool *has_ghosts );

/*! Get the data for a polyhedron list. */
typedef void ( *DTK_PolyhedronListDataFunction )(
    void *user_data, Coordinate *coordinates, LocalOrdinal *faces,
    unsigned *nodes_per_face, LocalOrdinal *cells, unsigned *faces_per_cell,
    int *face_orientation, bool *is_ghost_cell );

/*! Get the size parameters for building a cell list with a single topology. */
typedef void ( *DTK_CellListSizeFunction )(
    void *user_data, unsigned *space_dim, size_t *local_num_nodes,
    size_t *local_num_cells, unsigned *nodes_per_cell, bool *has_ghosts );

/*! Get the data for a single topology cell list. */
typedef void ( *DTK_CellListDataFunction )( void *user_data,
                                            Coordinate *coordinates,
                                            LocalOrdinal *cells,
                                            bool *is_ghost_cell,
                                            char *cell_topology );

/*! Get the size parameters for building a cell list with mixed topologies. */
typedef void ( *DTK_MixedTopologyCellListSizeFunction )(
    void *user_data, unsigned *space_dim, size_t *local_num_nodes,
    size_t *local_num_cells, size_t *total_nodes_per_cell, bool *has_ghosts );

/*! Get the data for a mixed topology cell list. */
typedef void ( *DTK_MixedTopologyCellListDataFunction )(
    void *user_data, Coordinate *coordinates, LocalOrdinal *cells,
    unsigned *cell_topology_ids, bool *is_ghost_cell, char **cell_topologies );

/*! Get the size parameters for a boundary. */
typedef void ( *DTK_BoundarySizeFunction )( void *user_data,
                                            const char *boundary_name,
                                            size_t *local_num_faces );

/*! Get the data for a boundary. */
typedef void ( *DTK_BoundaryDataFunction )( void *user_data,
                                            const char *boundary_name,
                                            LocalOrdinal *boundary_cells,
                                            unsigned *cell_faces_on_boundary );

/*! Get the size parameters for a degree-of-freedom id map with a single number
 * of dofs per object.
 */
typedef void ( *DTK_DOFMapSizeFunction )( void *user_data,
                                          size_t *local_num_dofs,
                                          size_t *local_num_objects,
                                          unsigned *dofs_per_object );

/*! Get the data for a degree-of-freedom id map with a single number of dofs
 * per object.
 */
typedef void ( *DTK_DOFMapDataFunction )( void *user_data,
                                          GlobalOrdinal *global_dof_ids,
                                          LocalOrdinal *object_dof_ids,
                                          char *discretization_type );

/*! Get the size parameters for a degree-of-freedom id map with each object
 * having a potentially different number of dofs (e.g. mixed topology cell
 * lists or polyhedron lists).
 */
typedef void ( *DTK_MixedTopologyDofMapSizeFunction )(
    void *user_data, size_t *local_num_dofs, size_t *local_num_objects,
    size_t *total_dofs_per_object );

/*! Get the data for a multiple object degree-of-freedom id map (e.g. mixed
 * topology cell lists or polyhedron lists).
 */
typedef void ( *DTK_MixedTopologyDofMapDataFunction )(
    void *user_data, GlobalOrdinal *global_dof_ids,
    LocalOrdinal *object_dof_ids, unsigned *dofs_per_object,
    char *discretization_type );

/*! Get the size parameters for a field. Field must be of size local_num_dofs
 * in the associated dof_id_map.
 */
typedef void ( *DTK_FieldSizeFunction )( void *user_data,
                                         const char *field_name,
                                         unsigned *field_dimension,
                                         size_t *local_num_dofs );

/*! Pull data from application into a field.  */
// FIXME: changed Scalar to double, which other do we need to provide?
typedef void ( *DTK_PullFieldDataFunction )( void *user_data,
                                             const char *field_name,
                                             double *field_dofs );

/* Push data from a field into the application. */
// FIXME: changed Scalar to double, which other do we need to provide?
typedef void ( *DTK_PushFieldDataFunction )( void *user_data,
                                             const char *field_name,
                                             const double *field_dofs );

/* Evaluate a field at a given set of points in a given set of objects. */
// FIXME: changed Scalar to double, which other do we need to provide?
typedef void ( *DTK_EvaluateFieldFunction )(
    void *user_data, const char *field_name,
    const Coordinate *evaluation_points, const LocalOrdinal *object_ids,
    double *values );

#ifdef __cplusplus
}
#endif

#endif // DTK_C_API_H
