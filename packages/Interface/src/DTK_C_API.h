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
#ifndef DTK_C_API_H
#define DTK_C_API_H

#include <DataTransferKit_config.hpp>

#include "DTK_CellTypes.h"
#include <DTK_Types.h>

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
extern const char *DTK_gitCommitHash();

/**
 * \defgroup c_interface_to_dtk_core Initialize/finalize DTK
 * @{
 */

/** \brief Initializes the DTK execution environment without any arguments.
 *
 *  This initializes Kokkos if it has not already been initialized.
 */
extern void DTK_initialize();

/** \brief Initialize DTK with command line arguments.
 *
 *  This initializes Kokkos if it has not already been initialized.
 *
 *  \param argc Pointer to the number of argument.
 *  \param argv Pointer to the argument vector.
 *
 *  This version of initialize() effectively calls <code>Kokkos::initialize(
 *  *argc, *argv )</code>.  Pointers to \c argc and \c argv arguments are passed
 *  in order to match MPI_Init's interface.
 */
extern void DTK_initializeCmd( int *argc, char ***argv );

/** \brief Indicates whether DTK has been initialized.
 *
 *  This function may be used to determine whether DTK has been initialized.
 */
extern bool DTK_isInitialized();

/** \brief Finalizes DTK.
 *
 *  This function terminates the DTK execution environment.  If DTK
 *  initialized Kokkos, this also finalizes Kokkos.  However, if Kokkos was
 *  initialized before DTK, then this function does NOT finalize Kokkos.
 */
extern void DTK_finalize();

/**@}*/

/**
 * \defgroup c_error_handling Error Handling
 * @{
 */

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
 * All DTK functions set \c errno error code upon completion. If a DTK
 * function fails, the code is nonzero. This function provides a way to get
 * the error message associated with the error code. If the error code is
 * unknown, DTK returns a string stating that.
 *
 * \param err Error number (typically, errno)
 * \return Returns corresponding error string. If the error code is 0 (success),
 *         return empty string.
 */
extern const char *DTK_error( int err );

/**@}*/

/**
 * \defgroup c_portability_enums Portability enumerations
 * @{
 */

/**
 *  \brief Memory space (where DTK memory is allocated)
 *
 *  A memory space defines the location in which DTK will allocate memory for
 *  the user to access. DTK callback functions allocate memory for reading
 *  inputs and writing outputs. Users declare the memory space in which they
 *  want their input and output arrays allocated via this enumeration. The
 *  following are valid values for the memory space enumeration:
 *
 *  DTK_HOST_SPACE: Memory will be allocated on the host in main CPU
 *  memory. Memory allocated in this space is not directly accessible on the
 *  GPU. If DTK maps are created for and executed on the GPU or other
 *  accelerators this memory will be explicitly copied to and from a memory
 *  space accessible by the GPU for computation.
 *
 *  DTK_CUDAUVM_SPACE: Memory will be allocated via CUDA unified virtual
 *  memory (UVM). This memory is automatically paged between host and device
 *  allowing users to access this memory both in standard host functions as
 *  well as in CUDA global and device functions. It should be noted that if
 *  this memory is accessed on the host in an offloading type of scenario then
 *  this memory will be paged between host and device when the user accesses
 *  the memory on the host and when DTK access the memory on the GPU and thus
 *  may incur a performance cost.
*/
typedef enum { DTK_HOST_SPACE, DTK_CUDAUVM_SPACE } DTK_MemorySpace;

/**
 *  \brief Execution space (where DTK functions execute)
 *
 *  An execution space defines where DTK mapping computations will occur on a
 *  compute node. Interpolation, projection, and other mathematical operations
 *  on user fields will execute using the programming model/runtime associated
 *  with the execution space parameter. The following are valid values for the
 *  execution space enumeration:
 *
 *  DTK_SERIAL: DTK kernels will execute in serial on a single CPU thread.
 *
 *  DTK_OPENMP: DTK kernels will execute in parallel on a number of OpenMP
 *  threads specified either via the environment variable OMP_NUM_THREADS or
 *  via specification to the Kokkos runtime in initialization via
 *  --kokkos-threads. If kokkos-specific runtime variables are used to
 *  specify thread counts, these should be passed at the time on DTK
 *  initialization using DTK_initializeCmd().
 *
 *  DTK_CUDA: DTK kernels will execute in parallel on an NVIDIA GPU using the
 *  CUDA runtime.
*/
typedef enum { DTK_SERIAL, DTK_OPENMP, DTK_CUDA } DTK_ExecutionSpace;

/**@}*/

/**
 * \defgroup c_interface_to_user_application User application data interface.
 * @{
 */

/**
 *  \brief DTK user application handle.
 *
 *  The user application handle represents an instance of the data access
 *  interface to a user application with user implementations of DTK call back
 *  functions are associated with each individual handle. As many handles may
 *  be created as desired with each representing its own unique instance.
 *
 *  \note the use of this handle in many interface functions below - all DTK
 *  functions needed access to user inputs and outputs will have user
 *  application handles as arguments.
 */
typedef struct _DTK_UserApplicationHandle *DTK_UserApplicationHandle;

/** \brief Create a DTK handle to a user application in a given memory space.
 *
 *  As many handles may be created as desired with each call to this function
 *  giving a new and unique handle. All data for user inputs and outputs
 *  accessed through function call backs registered with a given handle will
 *  be allocated in the memory space assocated with that handle. A call to
 *  this function should be associated with an equivalent call to
 *  DTK_destroyUserApplication when the handle's lifetime in the program is
 *  complete.
 *
 *  \note User application handles must be valid on all ranks of the
 *  communicator over which transfer operators are generated. In other words,
 *  this function must be called collectively on every rank in that
 *  communicator. In the case where the user's actual application does not
 *  exist on a given MPI rank (and therefore is represented by null data),
 *  this function must still be called. User implementations of call back
 *  functions for applications where the the actual application does not exist
 *  should simply return a size of zero for all application functions.
 *
 *  \param space Execution space for the callback functions that are to be
 *  registered using DTK_setUserFunction().
 *
 *  \return DTK_create returns a handle for the user application. All user
 *  inputs and outputs accessed through function call backs associated with
 *  this handle will be allocated in the given memory space.
 */
extern DTK_UserApplicationHandle
DTK_createUserApplication( DTK_MemorySpace space );

/** \brief Indicates whether a DTK handle to a user application is valid.
 *
 *  A handle is valid if it was created by DTK_createUserApplication() and has
 *  not yet been deleted by DTK_destroyUserApplication().
 *
 *  \param[in] handle The DTK user application handle to check.
 *
 *  \return true if the given user application handle is valid; false
 *  otherwise.
 */
extern bool DTK_isValidUserApplication( DTK_UserApplicationHandle handle );

/** \brief Destroy a DTK handle to a user application.
 *
 *  \param[in,out] handle User application handle. If this handle has already
 *  been destroyed or was not created with a call to
 *  DTK_createUserApplication() then this function does nothing.
 */
extern void DTK_destroyUserApplication( DTK_UserApplicationHandle handle );

/**@}*/

/**
 * \defgroup c_interface_to_map Interface to maps.
 * @{
 */

/** \brief DTK map handle.
 *
 *  The map handle represents a unique instance of a DTK transfer operator. A
 *  DTK map transfers data between a source (the application providing the
 *  data) and a target (the application receiving the data), each represented
 *  by their own user application handle providing access to their geometry,
 *  mesh, and field data. The type of map represented by the handle is defined
 *  via a set of input options and the execution space where map computations
 *  occur is defined at the time of construction via an execution space
 *  enumeration.
 *
 *  As many map instances may be created as needed and each instance may
 *  represent a different type of transfer operator. Once created, a map
 *  instance may be applied to transfer between the source and target as many
 *  times as needed as long as the source and target user application handles
 *  remain valid.
 */
typedef struct _DTK_MapHandle *DTK_MapHandle;

/** \brief Create a DTK handle to a transfer operator.
 *
 *  An options string is used to select the right map and parameters. This
 *  string is defined using a JSON syntax with key-value pairs. For example,
 *  specifying a nearest neighbor map would be achieved via:
 *  \code{.cpp}
 *      char *options = "{ \"Map Type\": \"Nearest Neighbor\" }";
 *  \endcode
 *  Some maps may have many options. These are separated via commas and may be
 *  strings or other types of plain-old-data. For example, consider making a
 *  map with an integer-valued option and a double-valued option:
 *  \code{.cpp}
 *      char *options = "{ \"Map Type\": \"Name Of Map\", "
 *                        "\"OptionFooInt\": 3, "
 *                        "\"OptionBarDouble\": 1.32 }";
 *  \endcode
 *
 *  \param space Execution space where the map will execute. Operations on
 *  user data for transfer operations will occur in this execution space. If
 *  the source or target applications reside in memory spaces that are not
 *  compatible with this execution space the data will be copied to and from a
 *  compatible memory space as needed.
 *
 *  \param[in] comm The MPI communicator over which to build the map.  Calls
 *  to both DTK_createMap() and DTK_applyMap() should be considered collective
 *  communications over this communicator. Note that this communicator must
 *  span all of the MPI ranks on which source and target data must be
 *  accessed. For example, if the source and target live on the same MPI
 *  communicator and therefore the same set of MPI ranks then that
 *  communicator should be the one passed to this function assuming that data
 *  for solution transfer will be accessed on all MPI ranks. If the source and
 *  target applications live on different MPI communicators composed of
 *  entirely different sets of MPI ranks then a new communicator that consists
 *  of all of the MPI ranks in both source and target communicators should be
 *  created and passed to this function. Cases will also arise in which the
 *  source or target application may not exist on some ranks of this
 *  communicator (e.g. the previously mentioned case of disjoint source and
 *  target communicators). In that case, user implementations of call back
 *  functions should just return sizes of zero during calls to allocation
 *  functions to indicate to DTK that there is no data from the user
 *  application on a given MPI rank.
 *
 *  \param[in] source Handle to the source application. This handle must be
 *  valid on all ranks in the communicator. Data will be pulled from this
 *  application and transferred to the target. Function call back
 *  implementations for the source should return zero sizes in allocation
 *  functions if the user's source application does not exist on the calling
 *  MPI rank.
 *
 *  \param[in,out] target Handle to the target application. Data will be
 *  transferred from the source and pushed to this application. This handle
 *  must be valid on all ranks in the communicator. Data will be pulled from
 *  this application and transferred to the target. Function call back
 *  implementations for the target should return zero sizes in allocation
 *  functions if the user's target application does not exist on the calling
 *  MPI rank.
 *
 *  \param[in] options Options string for building the map. See above for
 *  details on composing this option string.
 *
 *  \return DTK_create returns a handle for the map. This handle must be
 *  destroyed with DTK_destroyMap() when the lifetime of this map has ended in
 *  the program.
 */
extern DTK_MapHandle DTK_createMap( DTK_ExecutionSpace space, MPI_Comm comm,
                                    DTK_UserApplicationHandle source,
                                    DTK_UserApplicationHandle target,
                                    const char *options );

/** \brief Indicates whether a DTK handle to a map is valid.
 *
 *  A handle is valid if it was created by DTK_create() and has not yet been
 *  deleted by DTK_destroyUserApplication().
 *
 *  \param[in] handle The DTK map handle to check.
 *
 *  \return true if the given map handle is valid; false otherwise.
 */
extern bool DTK_isValidMap( DTK_MapHandle handle );

/** \brief Apply the DTK map to the given fields.
 *
 *  This function transfers the data from the source user application to the
 *  target user application. The fields transferred by this function are
 *  indicated by their given names. In practice, an application could
 *  implement their field function call backs to handle multiple fields,
 *  thereby allowing the same map instance to transfer many different fields.
 *
 *  \note This function call is a collective over the Map's communicator.
 *
 *  \note The source and target user application handles associated with the
 *  given map instance must still be valid - they cannot have been destroyed
 *  before this function is called.
 *
 *  \param[in] handle Map handle. This handle must be valid on all calling MPI
 *  ranks.
 *
 *  \param[in] source_field Name of the field in the source user
 *  application. This is the field name passed to either
 *  DTK_PullFieldDataFunction() or DTK_EvaluateFieldFunction() depending on
 *  the map type and user implementation of the source user application. DTK
 *  will read data from this field.
 *
 *  \param[in] target_field Name of the field in the target user
 *  application. This is the field name passed to DTK_PushFieldDataFunction()
 *  in the target user application. DTK will write data to this field.
 */
extern void DTK_applyMap( DTK_MapHandle handle, const char *source_field,
                          const char *target_field );

/** \brief Destroy a DTK handle to a map.
 *
 *  \param[in,out] handle map handle. If this handle has already been
 *  destroyed or was not created with a call to DTK_createMap() then this
 *  function does nothing.
 */
extern void DTK_destroyMap( DTK_MapHandle handle );

/**@}*/

/**
 * \defgroup c_function_registration User application data interface function
 * registration.
 * @{
 */

// clang-format off
// COMMENT: disabling clang-format because it keeps trying to put the comma on a
// separate new line.

/** \brief Enumeration passed as the \p type argument to DTK_setUserFunction()
 *  in order to indicate what callback function is being registered with the
 *  user application.
 *
 *  \note Callback functions are passed as pointers to functions that take no
 *  arguments and return nothing (<code>void(*)()</code>) so the value of the
 *  DTK_FunctionType enum is necessary to indicate which user function
 *  implementation is being registered with the user application interface.
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
// clang-format on

/** \brief Register a function as a callback.
 *
 *  This registers a custom user-implemented function as a callback for DTK to
 *  communicate with the user application. The registered function implements
 *  the interface specified by the function prototype (see prototypes below)
 *  associated with the given function type enumeration. When a function is
 *  registered, a user has an opportunity to also assign user data to that
 *  instance of the function via a <code>void*</code>. When the register user
 *  function is called, the pointer to user data is passed back to the
 *  function, allowing the user to use additional data as needed in the
 *  implementation or to customize the implementation for specific instances
 *  of their application.
 *
 *  \note Use the \p user_data pointer to your advantage when implementing
 *  user functions. Anything can be assigned to this pointer: a pointer
 *  directly to an instance of the actual user data, special data and
 *  functions written specifically for coupling, or other auxiliary data
 *  structures of the user's construction. Whatever you decide to put here
 *  will be passed back to you when the function is called - DTK will not
 *  modify this data whatsoever.
 *
 *  \param[in,out] handle User application handle. The function implementation
 *  will be registered with this particular instance.
 *  \param[in] type Type of callback function. The interface of the user
 *  function should correspond to the function prototype associated with this
 *  enumeration.
 *  \param[in] f Pointer to user defined callback function.
 *  \param[in] user_data Pointer to the user data that will be passed to the
 *             callback function when executing it.
 */
extern void DTK_setUserFunction( DTK_UserApplicationHandle handle,
                                 DTK_FunctionType type, void ( *f )(),
                                 void *user_data );

/**@}*/

/**
 * \defgroup c_function_prototypes User application data interface function
 * prototypes.
 * @{
 */

/** \brief Prototype function to get the size parameters for building a node
 *         list.
 *
 *  A node list is a collection of spatial points of a given dimension.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_NODE_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] space_dim Spatial dimension of the node coordinates.
 *  \param[out] local_num_nodes Number of nodes DTK will allocate memory for.
 */
typedef void ( *DTK_NodeListSizeFunction )( void *user_data,
                                            unsigned *space_dim,
                                            size_t *local_num_nodes );

/** \brief Prototype function to get the data for a node list.
 *
 *  A node is defined by its spatial coordinates.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_NODE_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] coordinates Node coordinates. Coordinates are blocked by
 *  dimension. For example, in 3 dimensions the x coordinates for all nodes
 *  are listed first followed by all of the y coordinates and then all of the
 *  z coordinates. A loop for this, for example, may look like:
 *  \code{.cpp}
 *      for ( int n = 0; n < local_num_nodes; ++n )
 *          for ( int d = 0; d < space_dim; ++d )
 *              coordinates[ d*local_num_nodes + n ] =
 *                  coordinate_of_node_n_in_dimension_d;
 *  \endcode
 */
typedef void ( *DTK_NodeListDataFunction )( void *user_data,
                                            Coordinate *coordinates );

/** \brief Prototype function to get the size parameters for building a bounding
 *  volume list.
 *
 *  A bounding volume list is a collection of axis-aligned Cartesian boxes in
 *  a given spatial dimension.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_BOUNDING_VOLUME_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *  \param[out] space_dim Spatial dimension.
 *  \param[out] local_num_volumes Number of volumes DTK will allocate memory
 *  for.
 */
typedef void ( *DTK_BoundingVolumeListSizeFunction )(
    void *user_data, unsigned *space_dim, size_t *local_num_volumes );

/** \brief Prototype function to get the data for a bounding volume list.
 *
 *  A bounding volume is defined by the low and high corner of the box
 *  (e.g. [x_min,y_min,z_min] and [x_max,y_max,z_max] in 3 dimensions).
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_BOUNDING_VOLUME_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] bounding_volumes Bounding volumes. This array specifies the
 *  coordinates of the low and high corners of each volume. The array is
 *  blocked by corner and the coordinates for each corner are blocked by
 *  dimension. The low corner comes first and the high corner comes second.
 *  A loop for this, for example, may look like:
 *  \code{.cpp}
 *      for ( int v = 0; v < local_num_volume; ++v )
 *          for ( int d = 0; d < space_dim; ++d )
 *          {
 *              bounding_volumes[ d*local_num_volumes + v ] =
 *                  low_corner_of_volume_v_in_dimension_d;
 *
 *              bounding_volumes[ (space_dim + d)*local_num_volumes + v ] =
 *                  high_corner_of_volume_v_in_dimension_d;
 *          }
 *  \endcode
 */
typedef void ( *DTK_BoundingVolumeListDataFunction )(
    void *user_data, Coordinate *bounding_volumes );

/** \brief Prototype function to get the size parameters for building a
 *  polyhedron list.
 *
 *  A polyehdron list is a collection of arbitrary linear polyhedra defined by
 *  a set of nodes and faces.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_POLYHEDRON_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] space_dim Spatial dimension.
 *
 *  \param[out] local_num_nodes Number of nodes DTK will allocate memory for.
 *
 *  \param[out] local_num_faces Number of faces DTK will allocate memory for.
 *
 *  \param[out] total_face_nodes Total number of nodes for all faces. This is
 *  equivalent to counting the number of nodes that construct each face and
 *  then summing this value over all faces.
 *
 *  \param[out] local_num_cells Number of cells DTK will allocate memory for.
 *
 *  \param[out] total_cell_faces Total number of faces for all cells. This is
 *  quivalent to counting the number of faces that construct each cell and
 *  then summing this value over all cells.
 */
typedef void ( *DTK_PolyhedronListSizeFunction )(
    void *user_data, unsigned *space_dim, size_t *local_num_nodes,
    size_t *local_num_faces, size_t *total_face_nodes, size_t *local_num_cells,
    size_t *total_cell_faces );

/** \brief Prototype function to get the data for a polyhedron list.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_POLYHEDRON_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] coordinates Node coordinates. This array is blocked by
 *  dimension in identical fashion to the NodeList coordinates.
 *
 *  \param[out] faces Connectivity list of faces. This array is defined as
 *  rank-1 but represents unstructured rank-2 data. It should be sized as
 *  (total sum of the number of nodes composing each face) or the sum of all
 *  elements in the following array, nodes_per_face, which indicates how many
 *  nodes are assigned to each face and how to index into this array. The
 *  input should be arranged as follows. Consider the \f$n^th\f$ node of face
 *  \f$i\f$ to be \f$f^i_n\f$ which is equal to the local index of the
 *  corresponding node in the coordinates array. Two faces, the first with 4
 *  nodes and the second with 3 would then be defined via this array as:
 *  \f$(f^1_1, f^1_2, f^1_3, f^1_4, f^2_1, f^2_2, f^2_3 )\f$ with the
 *  nodes_per_face array reading \f$(4, 3)\f$
 *
 *  \param[out] nodes_per_face Number of nodes per face. For every face, list
 *  how many nodes construct it. The sum of all local elements in this
 *  array should equal the total size of the faces array.
 *
 *  \param[out] cells Connectivity list of polyhedrons. This array is defined
 *  as rank-1 but represents unstructured rank-2 data. It should be sized as
 *  (total sum of the number of faces composing each polyhedron) or the sum of
 *  all elements in the array faces_per_cells, which indicates how many faces
 *  are assigned to each cell and how to index into this array. The input
 *  should be arranged as follows. Consider the \f$n^th\f$ face of cell
 *  \f$i\f$ to be \f$c^i_n\f$ which is equal to the local index of the
 *  corresponding face in the faces array. Two cells, the first with 5 faces
 *  and the second with 4 would then be defined via this array as: \f$(c^1_1,
 *  c^1_2, c^1_3, c^1_4, c^1_5, c^2_1, c^2_2, c^2_3, c^2_4 )\f$ with the
 *  faces_per_cell array reading \f$(5, 4)\f$.
 *
 *  \param[out] faces_per_cell Number of faces per cell. For every cell, list
 *  how many faces construct it. The sum of all local elements in this view
 *  should equal the total size of the cells view. This view is rank-1 and of
 *  length of the number of cells in the list.
 *
 *  \param[out] face_orientation Orientation of the faces.  Orientation of
 *  each face composing a cell indicating an outward or inward facing normal
 *  based on node ordering of the face and use of the right-hand rule. This
 *  view is defined as rank-1 but represents unstructured rank-2 data. This
 *  view is the same size as the cells view and is indexed in an identical
 *  matter. If the face for the given cell has a node ordering that returns a
 *  face normal that points into the cell via the right hand rule then a -1
 *  should be input. If the node ordering of the face produces a normal that
 *  points out from the cell a +1 should be input.
 */
typedef void ( *DTK_PolyhedronListDataFunction )(
    void *user_data, Coordinate *coordinates, LocalOrdinal *faces,
    unsigned *nodes_per_face, LocalOrdinal *cells, unsigned *faces_per_cell,
    int *face_orientation );

/** \brief Prototype function to get the size parameters for building a cell
 *  list.
 *
 *  Cells are objects from a topological zoo of cell types (e.g. hexahedron,
 *  triangle, etc.) and are defined by a cell type and a set of nodes ordered
 *  as prescribed by the cell type.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_CELL_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] space_dim Spatial dimension.
 *
 *  \param[out] local_num_nodes Number of nodes DTK will allocate memory for.
 *
 *  \param[out] local_num_cells Number of cells DTK will allocate memory for.
 *
 *  \param[out] total_cell_nodes Total number of nodes for all cells.
 */
typedef void ( *DTK_CellListSizeFunction )( void *user_data,
                                            unsigned *space_dim,
                                            size_t *local_num_nodes,
                                            size_t *local_num_cells,
                                            size_t *total_cell_nodes );

/** \brief Prototype function to get the data for a mixed topology cell list.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_CELL_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] coordinates Node coordinates. This array is blocked by
 *  dimension in identical fashion to the NodeList coordinates.
 *
 *  \param[out] cells List of cells. It represents a lists of cells with
 *  different topologies ordered in blocks. It should be sized as total sum of
 *  the number of nodes composing each cell. The input should be arranged as
 *  follows. Consider the \f$n^th\f$ node of cell \f$i\f$ to be \f$c^i_n\f$
 *  which is equal to the local index of the corresponding node in the nodes
 *  array. Two cells, the first with 5 nodes and the second with 4 would then
 *  be defined via this array as: \f$(c^1_1, c^1_2, c^1_3, c^1_4, c^1_5,
 *  c^2_1, c^2_2, c^2_3, c^2_4 )\f$ with the nodes_per_cell array reading
 *  \f$(5, 4)\f$. The number of nodes per cell is defined by the topology of
 *  the cell block given by the associated entry in block_topologies.
 *
 *  \param[out] cell_topologies Topologies of the cells. Give the cell
 *  topology type for each cell in the list.
 */
typedef void ( *DTK_CellListDataFunction )( void *user_data,
                                            Coordinate *coordinates,
                                            LocalOrdinal *cells,
                                            DTK_CellTopology *cell_topologies );

/** \brief Prototype function to get the size parameters for a boundary.
 *
 *  A boundary is a collection of cells (this includes faces of both
 *  polyhedrons and cells with fixed topologies) that coincide with a physical
 *  geometric boundary and the faces of those cells that are on the boundary.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_BOUNDARY_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] local_num_faces Number of faces owned by this process that are
 *  on the boundary.
 */
typedef void ( *DTK_BoundarySizeFunction )( void *user_data,
                                            size_t *local_num_faces );

/** \brief Prototype function to get the data for a boundary
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_BOUNDARY_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] boundary_cells Indices of the cells on the boundary. For every
 *  face on the boundary give the local id of the cell to which the face
 *  belongs. This array is of rank-1 and of length equal to the number of faces
 *  on the boundary. If the list does not have a boundary this array will be
 *  empty.
 *
 *  \param[out] cell_faces_on_boundary Indices of the faces within a given
 *  cell that is on the boundary. For every face on the boundary give the
 *  local id of the face relative to its parent cell. This is the local face
 *  id relative to the nodes as defined by the canonical cell topology. This
 *  array is of rank-1 and of length equal to the number of faces on the
 *  boundary. If the list does not have a boundary this array will be empty.
 */
typedef void ( *DTK_BoundaryDataFunction )( void *user_data,
                                            LocalOrdinal *boundary_cells,
                                            unsigned *cell_faces_on_boundary );

/** \brief Prototype function to get the size parameters for building an
 *  adjacency list.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_ADJACENCY_LIST_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] total_adjacencies Total number of adjacencies in the list.
 */
typedef void ( *DTK_AdjacencyListSizeFunction )( void *user_data,
                                                 size_t *total_adjacencies );

/** \brief Prototype function to get the data for an adjacency list.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_ADJACENCY_LIST_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] global_cell_ids The global ids of the local cells in the
 *  list.
 *
 *  \param[out] adjacent_global_cell_ids The global ids of the cells adjacent
 *  to the local cells in the list. These may live on another process.
 *
 *  \param[out] adjacencies_per_cell The number of adjacencies each local cell
 *  has. These serve as offsets into the adjacent_global_cell_ids array.
 */
typedef void ( *DTK_AdjacencyListDataFunction )(
    void *user_data, GlobalOrdinal *global_cell_ids,
    GlobalOrdinal *adjacent_global_cell_ids, unsigned *adjacencies_per_cell );

/** \brief Prototype function to get the size parameters for a
 *  degree-of-freedom id map with a single number of dofs per object.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_DOF_MAP_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] local_num_dofs Number of degrees of freedom owned by this
 *  process.
 *
 *  \param[out] local_num_objects Number of objects on this process.
 *
 *  \param[out] dofs_per_objects Degrees of freedom per object.
 */
typedef void ( *DTK_DOFMapSizeFunction )( void *user_data,
                                          size_t *local_num_dofs,
                                          size_t *local_num_objects,
                                          unsigned *dofs_per_object );

/** \brief Prototype function to get the size data for a degree-of-freedom id
 *  map with a single number of dofs per object.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_DOF_MAP_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] global_dof_ids Globally unique ids for DOFs on this process.
 *
 *  \param[out] object_dof_ids For every object of the given type in the object
 *  list give the local dof ids for that object. The local dof ids correspond
 *  to
 *
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
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] local_num_dofs Number of degrees of freedom owned by this
 *  process.
 *
 *  \param[out] local_num_objects Number of objects on this process.
 *
 *  \param[out] total_dofs_per_objects Total degrees of freedom per objects.
 */
typedef void ( *DTK_MixedTopologyDofMapSizeFunction )(
    void *user_data, size_t *local_num_dofs, size_t *local_num_objects,
    size_t *total_dofs_per_object );

/** \brief Prototype function to get the data for a multiple object
 *  degree-of-freedom id map (e.g. mixed topology cell lists or polyhedron
 *  lists).
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Pointer to custom user data.
 *
 *  \param[out] global_dof_ids Globally unique ids for DOFs on this process.
 *
 *  \param[out] object_dof_ids Local object IDs.
 *
 *  \param[out] dofs_per_object Degrees of freedom per object.
 *
 *  \param[out] discretization_type Type of discretization.
 */
typedef void ( *DTK_MixedTopologyDofMapDataFunction )(
    void *user_data, GlobalOrdinal *global_dof_ids,
    LocalOrdinal *object_dof_ids, unsigned *dofs_per_object,
    char *discretization_type );

/** \brief Prototype function to get the size parameters for a field.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_FIELD_SIZE_FUNCTION as the \p type argument.
 *
 *  Field must be of size local_num_dofs in the associated dof_id_map.
 *
 *  \param[in] user_data Custom user data.
 *
 *  \param[in] field_name Name of the field.
 *
 *  \param[in] field_dimension Dimension of the field (i.e. 1 for the pressure,
 *              or 3 for the velocity in 3-D)
 *
 *  \param[in] local_num_dofs Number of degrees of freedom owned by this
 *             process.
 */
typedef void ( *DTK_FieldSizeFunction )( void *user_data,
                                         const char *field_name,
                                         unsigned *field_dimension,
                                         size_t *local_num_dofs );

/** \brief Prototype function to pull data from the application into a field.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_PULL_FIELD_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Custom user data.
 *
 *  \param[in] field_name Name of the field to pull.
 *
 *  \param[out] field_dofs Degrees of freedom for that field.
 */
typedef void ( *DTK_PullFieldDataFunction )( void *user_data,
                                             const char *field_name,
                                             double *field_dofs );

/** \brief Prototype function to push data from a field into the application.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_PUSH_FIELD_DATA_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Custom user data.
 *
 *  \param[in] field_name Name of the field to push.
 *
 *  \param[out] field_dofs Degrees of freedom for that field.
 */
typedef void ( *DTK_PushFieldDataFunction )( void *user_data,
                                             const char *field_name,
                                             const double *field_dofs );

/** \brief Prototype function to evaluate a field at a given set of points in a
 *         given set of objects.
 *
 *  \note Register with a user application using DTK_setUserFunction() by
 *  passing DTK_EVALUATE_FIELD_FUNCTION as the \p type argument.
 *
 *  \param[in] user_data Custom user data.
 *
 *  \param[in] field_name Name of the field to evaluate.
 *
 *  \param[in] evaluate_points Coordinates of the points at which to evaluate
 *
 *             the field.
 *  \param[in] objects_ids ID of the cell/face with repect of which the
 *             coordinates are expressed.
 *
 *  \param[out] values Field values.
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
