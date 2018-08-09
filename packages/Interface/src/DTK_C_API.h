/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file
 * \brief C interface header.
 */
#ifndef DTK_C_API_H
#define DTK_C_API_H

#include <DataTransferKit_config.hpp>

#include <DTK_Types.h>
#include "DTK_CellTypes.h"

#include <mpi.h>

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
 *  Must be created using DTK_createUserApplication() to be a valid handle.
 *
 *  The handle essentially hides C++ implementation details from the user.
 *
 *  <!--
 *  Use incomplete types to differentiate between handles.
 *  We never define the incomplete structs.
 *  -->
 */
typedef struct _DTK_UserApplicationHandle *DTK_UserApplicationHandle;

/** \brief Memory space (where memory is allocated) */
typedef enum { DTK_HOST_SPACE, DTK_CUDAUVM_SPACE } DTK_MemorySpace;

/** \brief Execution space (where functions execute) */
typedef enum { DTK_SERIAL, DTK_OPENMP, DTK_CUDA } DTK_ExecutionSpace;

/** \brief Create a DTK handle to a user application.
 *
 *  \param space Execution space for the callback functions that are to be
 *  registered using DTK_setUserFunction().
 *
 *  \return DTK_create returns a handle for the user application.
 */
extern DTK_UserApplicationHandle
DTK_createUserApplication( DTK_MemorySpace space );

/** \brief Indicates whether a DTK handle to a user application is valid.
 *
 *  A handle is valid if it was created by DTK_create() and has not yet been
 *  deleted by DTK_destroy().
 *
 *  \param[in] handle The DTK user application handle to check.
 *
 *  \return true if the given user application handle is valid;  false
 *  otherwise.
 */
extern bool DTK_isValidUserApplication( DTK_UserApplicationHandle handle );

/** \brief Destroy a DTK handle to a user application.
 *
 *  \param[in,out] handle User application handle.
 */
extern void DTK_destroyUserApplication( DTK_UserApplicationHandle handle );

/**@}*/


/**
 * \defgroup c_interface_to_map Interface to maps.
 * @{
 */

/** \brief DTK map handle.
 *
 *  Must be created using DTK_createMap() to be a valid handle.
 *
 *  The handle essentially hides C++ implementation details from the user.
 *
 *  <!--
 *  Use incomplete types to differentiate between handles.
 *  We never define the incomplete structs.
 *  -->
 */
typedef struct _DTK_MapHandle *DTK_MapHandle;

/** \brief Create a DTK handle to a user appliction.
 *
 *  \param space Execution space where the map will execute.
 *
 *
 *  \param[in] comm The MPI communicator over which to build the map.
 *
 *  \param[in] source Handle to the source application.
 *
 *  \param[in,out] target Handle to the target application.
 *
 *  \return DTK_create returns a handle for the map.
 */
extern DTK_MapHandle DTK_createMap( DTK_ExecutionSpace space,
                                    MPI_Comm comm,
                                    DTK_UserApplicationHandle source,
                                    DTK_UserApplicationHandle target );

/** \brief Indicates whether a DTK handle to a map is valid.
 *
 *  A handle is valid if it was created by DTK_create() and has not yet been
 *  deleted by DTK_destroy().
 *
 *  \param[in] handle The DTK map handle to check.
 *
 *  \return true if the given map handle is valid; false otherwise.
 */
extern bool DTK_isValidMap( DTK_MapHandle handle );

/** \brief Apply the DTK map to the given fields.
 *
 *  \param[in] handle Map handle.
 *
 *  \param[in] source_field Name of the field in the source application.
 *
 *  \param[in] target_field Name of the field in the target application.
 */
extern void DTK_applyMap( DTK_MapHandle handle,
                          const char* source_field,
                          const char* target_field );

/** \brief Destroy a DTK handle to a map.
 *
 *  \param[in,out] handle map handle.
 */
extern void DTK_destroyMap( DTK_MapHandle handle );

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

/** \brief DTK error codes.
 */
typedef enum {
    DTK_SUCCESS = 0,
    DTK_INVALID_HANDLE = -1,
    DTK_UNINITIALIZED = -2,
    DTK_UNKNOWN = -99
} DTK_Error;

/** \brief Get DTK error message.
 *
 * All DTK functions set \c errno error code upon completion. If DTK function
 * fails, the code is nonzero. This function provides a way to get the the
 * error message associated with the error code. If the error code is unknown,
 * DTK returns a string stating that.
 *
 * \param err Error number (typically, errno)
 * \return Returns corresponding error string. If the error code is 0 (success),
 *         return empty string.
 */
extern const char *DTK_error( int err );

// clang-format off ////////////////////////////////////////////////////////////
// COMMENT: disabling clang-format because it keeps trying to put the comma on a
// separate new line.

/** \brief Passed as the \p type argument to DTK_setUserFunction() in order to
 *  indicate what callback function is being registered with the user application.
 *
 *  \note Callback functions are passed as pointers to functions that take no
 *  arguments and return nothing (<code>void(*)()</code>) so the value of the
 *  DTK_FunctionType enum is necessary to indicate what is being registered with
 *  the user application and how to cast the function pointer back to the
 *  appropriate signature.
 */
typedef enum {
    DTK_NODE_LIST_SIZE_FUNCTION /** See DTK_NodeListSizeFunction() */,
    DTK_NODE_LIST_DATA_FUNCTION /** See DTK_NodeListDataFunction() */,
    DTK_BOUNDING_VOLUME_LIST_SIZE_FUNCTION /** See DTK_BoundingVolumeListSizeFunction() */,
    DTK_BOUNDING_VOLUME_LIST_DATA_FUNCTION /** See DTK_BoundingVolumeListDataFunction() */,
    DTK_POLYHEDRON_LIST_SIZE_FUNCTION /** See DTK_PolyhedronListSizeFunction() */,
    DTK_POLYHEDRON_LIST_DATA_FUNCTION /** See DTK_PolyhedronListDataFunction() */,
    DTK_CELL_LIST_SIZE_FUNCTION /** See DTK_CellListSizeFunction() */,
    DTK_CELL_LIST_DATA_FUNCTION /** See DTK_CellListDataFunction() */,
    DTK_BOUNDARY_SIZE_FUNCTION /** See DTK_BoundarySizeFunction() */,
    DTK_BOUNDARY_DATA_FUNCTION /** See DTK_BoundaryDataFunction() */,
    DTK_ADJACENCY_LIST_SIZE_FUNCTION /** See DTK_AdjacencyListSizeFunction() */,
    DTK_ADJACENCY_LIST_DATA_FUNCTION /** See DTK_AdjacencyListDataFunction() */,
    DTK_DOF_MAP_SIZE_FUNCTION /** See DTK_DOFMapSizeFunction() */,
    DTK_DOF_MAP_DATA_FUNCTION /** See DTK_DOFMapDataFunction() */,
    DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION /** See DTK_MixedTopologyDofMapSizeFunction() */,
    DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION /** See DTK_MixedTopologyDofMapDataFunction() */,
    DTK_FIELD_SIZE_FUNCTION /** See DTK_FieldSizeFunction() */,
    DTK_PULL_FIELD_DATA_FUNCTION /** See DTK_PullFieldDataFunction() */,
    DTK_PUSH_FIELD_DATA_FUNCTION /** See DTK_PushFieldDataFunction() */,
    DTK_EVALUATE_FIELD_FUNCTION /** See DTK_EvaluateFieldFunction() */,
} DTK_FunctionType;
// clang-format on /////////////////////////////////////////////////////////////

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
extern void DTK_setUserFunction( DTK_UserApplicationHandle handle,
                                 DTK_FunctionType type, void ( *f )(),
                                 void *user_data );

/**
 * \defgroup c_interface_callbacks Prototype declaration of the callback
 * functions.
 * @{
 */

/** \brief Prototype function to get the size parameters for building a node
 *         list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_NODE_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] space_dim Spatial dimension.
 *  \param[out] local_num_nodes Number of nodes DTK will allocate memory for.
 */
typedef void ( *DTK_NodeListSizeFunction )( void *user_data,
                                            unsigned *space_dim,
                                            size_t *local_num_nodes );

/** \brief Prototype function to get the data for a node list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_NODE_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] coordinates Node coordinates.
 */
typedef void ( *DTK_NodeListDataFunction )( void *user_data,
                                            Coordinate *coordinates );

/** \brief Prototype function to get the size parameters for building a bounding
 *  volume list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_BOUNDING_VOLUME_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] space_dim Spatial dimension.
 *  \param[out] local_num_volumes Number of volumes DTK will allocate memory
 *  for.
 */
typedef void ( *DTK_BoundingVolumeListSizeFunction )( void *user_data,
                                                      unsigned *space_dim,
                                                      size_t *local_num_volumes );

/** \brief Prototype function to get the data for a bounding volume list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_BOUNDING_VOLUME_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] bounding_volumes Bounding volumes.
 */
typedef void ( *DTK_BoundingVolumeListDataFunction )(
    void *user_data, Coordinate *bounding_volumes );

/** \brief Prototype function to get the size parameters for building a
 *  polyhedron list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_POLYHEDRON_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] space_dim Spatial dimension.
 *  \param[out] local_num_nodes Number of nodes DTK will allocate memory for.
 *  \param[out] local_num_faces Number of faces DTK will allocate memory for.
 *  \param[out] total_face_nodes Total number of nodes for all faces.
 *  \param[out] local_num_cells Number of cells DTK will allocate memory for.
 *  \param[out] total_cell_faces Total number of faces for all cells.
 */
typedef void ( *DTK_PolyhedronListSizeFunction )(
    void *user_data, unsigned *space_dim, size_t *local_num_nodes,
    size_t *local_num_faces, size_t *total_face_nodes,
    size_t *local_num_cells, size_t *total_cell_faces );

/** \brief Prototype function to get the data for a polyhedron list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_POLYHEDRON_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] coordinates Node coordinates.
 *  \param[out] faces Connectivity list of faces.
 *  \param[out] nodes_per_face Number of nodes per face.
 *  \param[out] cells Connectivity list of polyhedrons.
 *  \param[out] faces_per_cell Number of faces per cell.
 *  \param[out] face_orientation Orientation of the faces.
 */
typedef void ( *DTK_PolyhedronListDataFunction )(
    void *user_data, Coordinate *coordinates, LocalOrdinal *faces,
    unsigned *nodes_per_face, LocalOrdinal *cells, unsigned *faces_per_cell,
    int *face_orientation );

/** \brief Prototype function to get the size parameters for building a cell
 *  list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_CELL_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] space_dim Spatial dimension.
 *  \param[out] local_num_nodes Number of nodes DTK will allocate memory for.
 *  \param[out] local_num_cells Number of cells DTK will allocate memory for.
 *  \param[out] total_cell_nodes Total number of nodes for all cells.
 */
typedef void ( *DTK_CellListSizeFunction )(
    void *user_data, unsigned *space_dim, size_t *local_num_nodes,
    size_t *local_num_cells, size_t *total_cell_nodes );

/** \brief Prototype function to get the data for a mixed topology cell list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_CELL_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] coordinates Node coordinates.
 *  \param[out] cells List of cells.
 *  \param[out] cell_topologies Topologies of the cells.
 */
typedef void ( *DTK_CellListDataFunction )(
    void *user_data, Coordinate *coordinates, LocalOrdinal *cells,
    DTK_CellTopology *cell_topologies );

/** \brief Prototype function to get the size parameters for a boundary
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_BOUNDARY_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] local_num_faces Number of faces owned by this process.
 */
typedef void ( *DTK_BoundarySizeFunction )( void *user_data,
                                            size_t *local_num_faces );

/** \brief Prototype function to get the data for a boundary
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_BOUNDARY_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] boundary_cells Indices of the cells on the boundary.
 *  \param[out] cell_faces_on_boundary Indices of the faces within a given cell
 *  that is on the boundary.
 */
typedef void ( *DTK_BoundaryDataFunction )( void *user_data,
                                            LocalOrdinal *boundary_cells,
                                            unsigned *cell_faces_on_boundary );

/** \brief Prototype function to get the size parameters for building an
 *  adjacency list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_ADJACENCY_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] total_adjacencies Total number of adjacencies in the list.
 */
typedef void ( *DTK_AdjacencyListSizeFunction )(
    void *user_data, size_t *total_adjacencies );

/** \brief Prototype function to get the data for an adjacency list.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_ADJACENCY_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] global_cell_ids The global ids of the local cells in the
 *  list.
 *  \param[out] adjacent_global_cell_ids The global ids of the cells adjacent
 *  to the local cells in the list. These may live on another process.
 *  \param[out] adjacencies_per_cell The number of adjacencies each local cell
 *  has. These serve as offsets into the adjacent_global_cell_ids array.
 */
typedef void ( *DTK_AdjacencyListDataFunction )(
    void *user_data, GlobalOrdinal* global_cell_ids,
    GlobalOrdinal* adjacent_global_cell_ids,
    unsigned* adjacencies_per_cell );

/** \brief Prototype function to get the size parameters for a
 *  degree-of-freedom id map with a single number of dofs per object.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_DOF_MAP_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] local_num_dofs Number of degrees of freedom owned by this
 *  process.
 *  \param[out] local_num_objects Number of objects on this process.
 *  \param[out] dofs_per_objects Degrees of freedom per object.
 */
typedef void ( *DTK_DOFMapSizeFunction )( void *user_data,
                                          size_t *local_num_dofs,
                                          size_t *local_num_objects,
                                          unsigned *dofs_per_object );

/** \brief Prototype function to get the size data for a degree-of-freedom id
 *  map with a single number of dofs per object.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_DOF_MAP_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] global_dof_ids Globally unique ids for DOFs on this process.
 *  \param[out] object_dof_ids For every object of the given type in the object
 *  list give the local dof ids for that object. The local dof ids correspond to
 *  the index of the entry in the global dof id view.
 *  \param[out] discretization_type Type of discretization.
 */
typedef void ( *DTK_DOFMapDataFunction )( void *user_data,
                                          GlobalOrdinal *global_dof_ids,
                                          LocalOrdinal *object_dof_ids,
                                          char *discretization_type );

/** \brief Prototype function to get the size parameters for a
 *  degree-of-freedom id map with each object having a potentially different
 *  number of dofs (e.g. mixed topology cell lists or polyhedron lists).
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] local_num_dofs Number of degrees of freedom owned by this
 *  process.
 *  \param[out] local_num_objects Number of objects on this process.
 *  \param[out] total_dofs_per_objects Total degrees of freedom per objects.
 */
typedef void ( *DTK_MixedTopologyDofMapSizeFunction )(
    void *user_data, size_t *local_num_dofs, size_t *local_num_objects,
    size_t *total_dofs_per_object );

/** \brief Prototype function to get the data for a multiple object
 *  degree-of-freedom id map (e.g. mixed topology cell lists or polyhedron
 *  lists).
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] global_dof_ids Globally unique ids for DOFs on this process.
 *  \param[out] object_dof_ids Local object IDs.
 *  \param[out] dofs_per_object Degrees of freedom per object.
 *  \param[out] discretization_type Type of discretization.
 */
typedef void ( *DTK_MixedTopologyDofMapDataFunction )(
    void *user_data, GlobalOrdinal *global_dof_ids,
    LocalOrdinal *object_dof_ids, unsigned *dofs_per_object,
    char *discretization_type );

/** \brief Prototype function to get the size parameters for a field.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_FIELD_SIZE_FUNCTION as the \p type argument.
 *
 *  Field must be of size local_num_dofs in the associated dof_id_map.
 *
 *  \param[in] user_data Custom user data.
 *  \param[in] field_name Name of the field.
 *  \param[in] field_dimension Dimension of the field (i.e. 1 for the pressure,
 *              or 3 for the velocity in 3-D)
 *  \param[in] local_num_dofs Number of degrees of freedom owned by this
 *             process.
 */
typedef void ( *DTK_FieldSizeFunction )( void *user_data,
                                         const char *field_name,
                                         unsigned *field_dimension,
                                         size_t *local_num_dofs );

/** \brief Prototype function to pull data from the application into a field.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_PULL_FIELD_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Custom user data.
 *  \param[in] field_name Name of the field to pull.
 *  \param[out] field_dofs Degrees of freedom for that field.
 *  <!--
 *  FIXME: changed Scalar to double, which other do we need to provide?
 *  -->
 */
typedef void ( *DTK_PullFieldDataFunction )( void *user_data,
                                             const char *field_name,
                                             double *field_dofs );

/** \brief Prototype function to push data from a field into the application.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_PUSH_FIELD_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Custom user data.
 *  \param[in] field_name Name of the field to push.
 *  \param[out] field_dofs Degrees of freedom for that field.
 *  <!--
 *  FIXME: changed Scalar to double, which other do we need to provide?
 *  -->
 */
typedef void ( *DTK_PushFieldDataFunction )( void *user_data,
                                             const char *field_name,
                                             const double *field_dofs );

/** \brief Prototype function to evaluate a field at a given set of points in a
 *         given set of objects.
 *
 *  Register with a user application using DTK_setUserFunction() by passing
 *  DTK_EVALUATE_FIELD_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Custom user data.
 *  \param[in] field_name Name of the field to evaluate.
 *  \param[in] evaluate_points Coordinates of the points at which to evaluate
 *             the field.
 *  \param[in] objects_ids ID of the cell/face with repect of which the
 *             coordinates are expressed.
 *  \param[out] values Field values.
 *  <!--
 *  FIXME: changed Scalar to double, which other do we need to provide?
 *  -->
 */
typedef void ( *DTK_EvaluateFieldFunction )(
    void *user_data, const char *field_name,
    const Coordinate *evaluation_points, const LocalOrdinal *object_ids,
    double *values );

/**@}*/

#ifdef __cplusplus
}
#endif

#endif // DTK_C_API_H
