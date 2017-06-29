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

/**
 * \defgroup c_interface_runtime_api C runtime API
 * @{
 */

/** \brief Get the current version of DTK.
 *
 *  \return Returns a string containing the version number for DTK.
 */
extern const char *DTK_version();

/** \brief Get the current repository hash.
 *
 *  \note If the source code is not under revision control (e.g. downloaded as
 *  a tarball), this functions returns an error string indicating that it is not
 *  a GIT repository.
 *
 *  \return Returns a string containing the revision number.
 */
extern const char *DTK_git_commit_hash();

/**@}*/

/**
 * \defgroup c_interface_to_user_application Interface to user app.
 * @{
 */

/** \brief DTK user application handle.
 *
 *  Must be created using DTK_create() to be a valid handle.
 *
 *  The handle essentially hides C++ implementation details from the user.
 *
 *  <!--
 *  Use incomplete types to differentiate between handles.
 *  We never define the incomplete structs.
 *  -->
 */
typedef struct _DTK_UserApplicationHandle *DTK_UserApplicationHandle;

/** \brief Execution space (where functions execute) */
typedef enum { DTK_SERIAL, DTK_OPENMP, DTK_CUDA } DTK_ExecutionSpace;

/** \brief Create a DTK handle.
 *
 *  \param space Execution space for the callback functions that are to be
 *  registered using DTK_set_function().
 *
 *  \return DTK_create returns a handle for the user application.
 */
extern DTK_UserApplicationHandle DTK_create( DTK_ExecutionSpace space );

/** \brief Indicates whether a DTK handle is valid.
 *
 *  A handle is valid if it was created by DTK_create() and has not yet been
 *  deleted by DTK_destroy().
 *
 *  \param[in] handle The DTK user application handle to check.
 *
 *  \return true if the given user application handle is valid;  false
 *  otherwise.
 */
extern bool DTK_is_valid( DTK_UserApplicationHandle handle );

/** \brief Destroy a DTK handle.
 *
 *  \param[in,out] handle User application handle.
 */
extern void DTK_destroy( DTK_UserApplicationHandle handle );

/**@}*/

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

/** \brief Passed as the \p type argument to DTK_set_function() in order to
 *  indicate what callback function is being registered with the user application.
 *
 *  \note Callback functions are passed as pointers to functions that take no
 *  arguments and return nothing (<code>void(*)()</code>) so the value of the
 *  DTK_FunctionType enum is necessary to indicate what is being registered with
 *  the user application and how to cast the function pointer back to the
 *  appropriate signature.
 */
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

/** \brief Register a function as a callback.
 *
 *  This registers a custom function as a callback for DTK to communicate with
 *  the user application.
 *
 *  \param[in,out] handle User application handle.
 *  \param[in] type Type of callback function.
 *  \param[in] f Pointer to user defined callback function.
 *  \param[in] user_data Pointer to the user data that will be passed to the
 *             callback function when executing it.
 */
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
