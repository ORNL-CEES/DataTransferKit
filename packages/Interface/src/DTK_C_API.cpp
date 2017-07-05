/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include "DTK_C_API.hpp"
#include "DTK_Core.hpp"

#include "DTK_Version.hpp"

#include <cerrno>
#include <set>

namespace DataTransferKit
{

struct DTK_FunctionWrap
{
    DTK_FunctionWrap( void ( *f )(), void *data )
        : _f( f )
        , _data( data )
    {
    }
    void ( *_f )();
    void *_data;
};

// We store the reinterpret_cast versions of pointers
static std::set<void *> valid_handles;

template <typename Function>
std::pair<Function, void *> get_function( std::shared_ptr<void> user_data )
{
    auto u = std::static_pointer_cast<DTK_FunctionWrap>( user_data );

    return std::make_pair( Function( u->_f ), u->_data );
}

void NodeListSizeFunctionWrapper( std::shared_ptr<void> user_data,
                                  unsigned &space_dim, size_t &local_num_nodes,
                                  bool &has_ghosts )
{
    auto u = get_function<DTK_NodeListSizeFunction>( user_data );
    u.first( u.second, &space_dim, &local_num_nodes, &has_ghosts );
}

void NodeListDataFunctionWrapper( std::shared_ptr<void> user_data,
                                  View<Coordinate> coordinates,
                                  View<bool> is_ghost_node )
{
    auto u = get_function<DTK_NodeListDataFunction>( user_data );
    u.first( u.second, coordinates.data(), is_ghost_node.data() );
}

void BoundingVolumeListSizeFunctionWrapper( std::shared_ptr<void> user_data,
                                            unsigned &space_dim,
                                            size_t &local_num_volumes,
                                            bool &has_ghosts )
{
    auto u = get_function<DTK_BoundingVolumeListSizeFunction>( user_data );
    u.first( u.second, &space_dim, &local_num_volumes, &has_ghosts );
}

void BoundingVolumeListDataFunctionWrapper( std::shared_ptr<void> user_data,
                                            View<Coordinate> bounding_volumes,
                                            View<bool> is_ghost_volume )
{
    auto u = get_function<DTK_BoundingVolumeListDataFunction>( user_data );
    u.first( u.second, bounding_volumes.data(), is_ghost_volume.data() );
}

void PolyhedronListSizeFunctionWrapper(
    std::shared_ptr<void> user_data, unsigned &space_dim,
    size_t &local_num_nodes, size_t &local_num_faces,
    size_t &total_nodes_per_face, size_t &local_num_cells,
    size_t &total_faces_per_cell, bool &has_ghosts )
{
    auto u = get_function<DTK_PolyhedronListSizeFunction>( user_data );
    u.first( u.second, &space_dim, &local_num_nodes, &local_num_faces,
             &total_nodes_per_face, &local_num_cells, &total_faces_per_cell,
             &has_ghosts );
}

void PolyhedronListDataFunctionWrapper(
    std::shared_ptr<void> user_data, View<Coordinate> coordinates,
    View<LocalOrdinal> faces, View<unsigned> nodes_per_face,
    View<LocalOrdinal> cells, View<unsigned> faces_per_cell,
    View<int> face_orientation, View<bool> is_ghost_cell )
{
    auto u = get_function<DTK_PolyhedronListDataFunction>( user_data );
    u.first( u.second, coordinates.data(), faces.data(), nodes_per_face.data(),
             cells.data(), faces_per_cell.data(), face_orientation.data(),
             is_ghost_cell.data() );
}

void CellListSizeFunctioWrapper( std::shared_ptr<void> user_data,
                                 unsigned &space_dim, size_t &local_num_nodes,
                                 size_t &local_num_cells,
                                 unsigned &nodes_per_cell, bool &has_ghosts )
{
    auto u = get_function<DTK_CellListSizeFunction>( user_data );
    u.first( u.second, &space_dim, &local_num_nodes, &local_num_cells,
             &nodes_per_cell, &has_ghosts );
}

void CellListDataFunctionWrapper( std::shared_ptr<void> user_data,
                                  View<Coordinate> coordinates,
                                  View<LocalOrdinal> cells,
                                  View<bool> is_ghost_cell,
                                  std::string &cell_topology )
{
    const int max_string_size = 255;
    std::vector<char> c_cell_topology( max_string_size );

    auto u = get_function<DTK_CellListDataFunction>( user_data );
    u.first( u.second, coordinates.data(), cells.data(), is_ghost_cell.data(),
             c_cell_topology.data() );

    cell_topology.assign( c_cell_topology.data() );
}

void MixedTopologyCellListSizeFunctionWrapper( std::shared_ptr<void> user_data,
                                               unsigned &space_dim,
                                               size_t &local_num_nodes,
                                               size_t &local_num_cells,
                                               size_t &total_nodes_per_cell,
                                               bool &has_ghosts )
{
    auto u = get_function<DTK_MixedTopologyCellListSizeFunction>( user_data );
    u.first( u.second, &space_dim, &local_num_nodes, &local_num_cells,
             &total_nodes_per_cell, &has_ghosts );
}

void MixedTopologyCellListDataFunctionWrapper(
    std::shared_ptr<void> user_data, View<Coordinate> coordinates,
    View<LocalOrdinal> cells, View<unsigned> cell_topology_ids,
    View<bool> is_ghost_cell, std::vector<std::string> &cell_topologies )
{
    const int max_num_topologies = 100;
    const int max_string_size = 255;

    std::vector<char *> c_cell_topologies( max_num_topologies );
    for ( int i = 0; i < max_num_topologies; i++ )
    {
        c_cell_topologies[i] = new char[max_string_size];
        c_cell_topologies[i][0] = 0;
    }

    auto u = get_function<DTK_MixedTopologyCellListDataFunction>( user_data );
    u.first( u.second, coordinates.data(), cells.data(),
             cell_topology_ids.data(), is_ghost_cell.data(),
             c_cell_topologies.data() );

    // count
    int cnt = 0;
    for ( ; cnt < max_num_topologies; cnt++ )
        if ( c_cell_topologies[cnt][0] == 0 )
            break;
    // assign
    cell_topologies.resize( cnt );
    for ( int i = 0; i < cnt; i++ )
    {
        cell_topologies[i].assign( c_cell_topologies[i] );
        delete[] c_cell_topologies[i];
    }
}

void BoundarySizeFunctionWrapper( std::shared_ptr<void> user_data,
                                  const std::string &boundary_name,
                                  size_t &local_num_faces )
{
    auto u = get_function<DTK_BoundarySizeFunction>( user_data );
    u.first( u.second, boundary_name.c_str(), &local_num_faces );
}

void BoundaryDataFunctionWrapper( std::shared_ptr<void> user_data,
                                  const std::string &boundary_name,
                                  View<LocalOrdinal> boundary_cells,
                                  View<unsigned> cell_faces_on_boundary )
{
    auto u = get_function<DTK_BoundaryDataFunction>( user_data );
    u.first( u.second, boundary_name.c_str(), boundary_cells.data(),
             cell_faces_on_boundary.data() );
}

void DOFMapSizeFunctionWrapper( std::shared_ptr<void> user_data,
                                size_t &local_num_dofs,
                                size_t &local_num_objects,
                                unsigned &dofs_per_object )
{
    auto u = get_function<DTK_DOFMapSizeFunction>( user_data );
    u.first( u.second, &local_num_dofs, &local_num_objects, &dofs_per_object );
}

void DOFMapDataFunctionWrapper( std::shared_ptr<void> user_data,
                                View<GlobalOrdinal> global_dof_ids,
                                View<LocalOrdinal> object_dof_ids,
                                std::string &discretization_type )
{
    const int max_string_size = 255;
    std::vector<char> c_discretization_type( max_string_size );

    auto u = get_function<DTK_DOFMapDataFunction>( user_data );
    u.first( u.second, global_dof_ids.data(), object_dof_ids.data(),
             c_discretization_type.data() );

    discretization_type.assign( c_discretization_type.data() );
}

void MixedTopologyDOFMapSizeFunctionWrapper( std::shared_ptr<void> user_data,
                                             size_t &local_num_dofs,
                                             size_t &local_num_objects,
                                             size_t &total_dofs_per_object )
{
    auto u = get_function<DTK_MixedTopologyDofMapSizeFunction>( user_data );
    u.first( u.second, &local_num_dofs, &local_num_objects,
             &total_dofs_per_object );
}

void MixedTopologyDOFMapDataFunctionWrapper( std::shared_ptr<void> user_data,
                                             View<GlobalOrdinal> global_dof_ids,
                                             View<LocalOrdinal> object_dof_ids,
                                             View<unsigned> dofs_per_object,
                                             std::string &discretization_type )
{
    const int max_string_size = 255;
    std::vector<char> c_discretization_type( max_string_size );

    auto u = get_function<DTK_MixedTopologyDofMapDataFunction>( user_data );
    u.first( u.second, global_dof_ids.data(), object_dof_ids.data(),
             dofs_per_object.data(), c_discretization_type.data() );

    discretization_type.assign( c_discretization_type.data() );
}

void FieldSizeFunctionWrapper( std::shared_ptr<void> user_data,
                               const std::string &field_name,
                               unsigned &field_dimension,
                               size_t &local_num_dofs )
{
    auto u = get_function<DTK_FieldSizeFunction>( user_data );
    u.first( u.second, field_name.c_str(), &field_dimension, &local_num_dofs );
}

template <class Scalar>
void PullFieldDataFunctionWrapper( std::shared_ptr<void> user_data,
                                   const std::string &field_name,
                                   View<Scalar> field_dofs )
{
    throw DataTransferKitException( "Not implemented" );
}
template <>
void PullFieldDataFunctionWrapper<double>( std::shared_ptr<void> user_data,
                                           const std::string &field_name,
                                           View<double> field_dofs )
{
    auto u = get_function<DTK_PullFieldDataFunction>( user_data );
    u.first( u.second, field_name.c_str(), field_dofs.data() );
}

template <class Scalar>
void PushFieldDataFunctionWrapper( std::shared_ptr<void> user_data,
                                   const std::string &field_name,
                                   const View<Scalar> field_dofs )
{
    throw DataTransferKitException( "Not implemented" );
}
template <>
void PushFieldDataFunctionWrapper<double>( std::shared_ptr<void> user_data,
                                           const std::string &field_name,
                                           const View<double> field_dofs )
{
    auto u = get_function<DTK_PushFieldDataFunction>( user_data );
    u.first( u.second, field_name.c_str(), field_dofs.data() );
}

template <class Scalar>
void EvaluateFieldFunctionWrapper( std::shared_ptr<void> user_data,
                                   const std::string &field_name,
                                   const View<Coordinate> evaluation_points,
                                   const View<LocalOrdinal> object_ids,
                                   View<Scalar> values )
{
    throw DataTransferKitException( "Not implemented" );
}
template <>
void EvaluateFieldFunctionWrapper<double>(
    std::shared_ptr<void> user_data, const std::string &field_name,
    const View<Coordinate> evaluation_points,
    const View<LocalOrdinal> object_ids, View<double> values )
{
    auto u = get_function<DTK_EvaluateFieldFunction>( user_data );
    u.first( u.second, field_name.c_str(), evaluation_points.data(),
             object_ids.data(), values.data() );
}

} // namespace DataTransferKit

extern "C" {

typedef enum {
    DTK_SUCCESS = 0,
    DTK_INVALID_HANDLE = -1,
    DTK_UNINITIALIZED = -2,
    DTK_UNKNOWN = -99
} DTK_Error;

const char *DTK_version()
{
    errno = DTK_SUCCESS;
    static std::string version_string = DataTransferKit::version().c_str();
    return version_string.c_str();
}

const char *DTK_git_commit_hash()
{
    errno = DTK_SUCCESS;
    static std::string hash_string = DataTransferKit::git_commit_hash().c_str();
    return hash_string.c_str();
}

DTK_UserApplicationHandle DTK_create( DTK_ExecutionSpace space )
{
    errno = DTK_SUCCESS;
    if ( !DTK_is_initialized() )
    {
        errno = DTK_UNINITIALIZED;
        return nullptr;
    }

    auto handle = reinterpret_cast<DTK_UserApplicationHandle>(
        new DataTransferKit::DTK_Registry( space ) );
    DataTransferKit::valid_handles.insert( handle );

    return handle;
}

bool DTK_is_valid( DTK_UserApplicationHandle handle )
{
    errno = DTK_SUCCESS;
    return DataTransferKit::valid_handles.count( handle );
}

void DTK_destroy( DTK_UserApplicationHandle handle )
{
    errno = DTK_SUCCESS;
    if ( DataTransferKit::valid_handles.count( handle ) )
    {
        auto dtk = reinterpret_cast<DataTransferKit::DTK_Registry *>( handle );
        // nullptr is definitely not a valid handle, so no need to check
        delete dtk;
        // use handle instead of dtk as reinterpret_cast may change pointers
        DataTransferKit::valid_handles.erase( handle );
    }
}

void DTK_initialize()
{
    errno = DTK_SUCCESS;
    DataTransferKit::initialize();
}

void DTK_initialize_cmd( int *argc, char ***argv )
{
    errno = DTK_SUCCESS;
    DataTransferKit::initialize( *argc, *argv );
}

bool DTK_is_initialized()
{
    errno = DTK_SUCCESS;
    return DataTransferKit::isInitialized();
}

void DTK_finalize()
{
    errno = DTK_SUCCESS;
    DataTransferKit::finalize();
}

void DTK_set_function( DTK_UserApplicationHandle handle, DTK_FunctionType type,
                       void ( *f )(), void *user_data )
{
    errno = DTK_SUCCESS;

    using namespace DataTransferKit;

    if ( !DTK_is_valid( handle ) )
    {
        errno = DTK_INVALID_HANDLE;
        return;
    }

    try
    {
        auto dtk = reinterpret_cast<DTK_Registry *>( handle );
        auto data = std::make_shared<DTK_FunctionWrap>( f, user_data );

        switch ( type )
        {
        case DTK_NODE_LIST_SIZE_FUNCTION:
            dtk->_registry->setNodeListSizeFunction(
                NodeListSizeFunctionWrapper, data );
            break;
        case DTK_NODE_LIST_DATA_FUNCTION:
            dtk->_registry->setNodeListDataFunction(
                NodeListDataFunctionWrapper, data );
            break;
        case DTK_BOUNDING_VOLUME_LIST_SIZE_FUNCTION:
            dtk->_registry->setBoundingVolumeListSizeFunction(
                BoundingVolumeListSizeFunctionWrapper, data );
            break;
        case DTK_BOUNDING_VOLUME_LIST_DATA_FUNCTION:
            dtk->_registry->setBoundingVolumeListDataFunction(
                BoundingVolumeListDataFunctionWrapper, data );
            break;
        case DTK_POLYHEDRON_LIST_SIZE_FUNCTION:
            dtk->_registry->setPolyhedronListSizeFunction(
                PolyhedronListSizeFunctionWrapper, data );
            break;
        case DTK_POLYHEDRON_LIST_DATA_FUNCTION:
            dtk->_registry->setPolyhedronListDataFunction(
                PolyhedronListDataFunctionWrapper, data );
            break;
        case DTK_CELL_LIST_SIZE_FUNCTION:
            dtk->_registry->setCellListSizeFunction( CellListSizeFunctioWrapper,
                                                     data );
            break;
        case DTK_CELL_LIST_DATA_FUNCTION:
            dtk->_registry->setCellListDataFunction(
                CellListDataFunctionWrapper, data );
            break;
        case DTK_MIXED_TOPOLOGY_CELL_LIST_SIZE_FUNCTION:
            dtk->_registry->setMixedTopologyCellListSizeFunction(
                MixedTopologyCellListSizeFunctionWrapper, data );
            break;
        case DTK_MIXED_TOPOLOGY_CELL_LIST_DATA_FUNCTION:
            dtk->_registry->setMixedTopologyCellListDataFunction(
                MixedTopologyCellListDataFunctionWrapper, data );
            break;
        case DTK_BOUNDARY_SIZE_FUNCTION:
            dtk->_registry->setBoundarySizeFunction(
                BoundarySizeFunctionWrapper, data );
            break;
        case DTK_BOUNDARY_DATA_FUNCTION:
            dtk->_registry->setBoundaryDataFunction(
                BoundaryDataFunctionWrapper, data );
            break;
        case DTK_DOF_MAP_SIZE_FUNCTION:
            dtk->_registry->setDOFMapSizeFunction( DOFMapSizeFunctionWrapper,
                                                   data );
            break;
        case DTK_DOF_MAP_DATA_FUNCTION:
            dtk->_registry->setDOFMapDataFunction( DOFMapDataFunctionWrapper,
                                                   data );
            break;
        case DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION:
            dtk->_registry->setMixedTopologyDOFMapSizeFunction(
                MixedTopologyDOFMapSizeFunctionWrapper, data );
            break;
        case DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION:
            dtk->_registry->setMixedTopologyDOFMapDataFunction(
                MixedTopologyDOFMapDataFunctionWrapper, data );
            break;
        case DTK_FIELD_SIZE_FUNCTION:
            dtk->_registry->setFieldSizeFunction( FieldSizeFunctionWrapper,
                                                  data );
            break;
        case DTK_PULL_FIELD_DATA_FUNCTION:
            dtk->_registry->setPullFieldDataFunction(
                PullFieldDataFunctionWrapper<double>, data );
            break;
        case DTK_PUSH_FIELD_DATA_FUNCTION:
            dtk->_registry->setPushFieldDataFunction(
                PushFieldDataFunctionWrapper<double>, data );
            break;
        case DTK_EVALUATE_FIELD_FUNCTION:
            dtk->_registry->setEvaluateFieldFunction(
                EvaluateFieldFunctionWrapper<double>, data );
        }
    }
    catch ( ... )
    {
        errno = DTK_UNKNOWN;
    }
}

const char *DTK_error( int err )
{
    errno = DTK_SUCCESS;
    switch ( err )
    {
    case DTK_SUCCESS:
        return "";
    case DTK_INVALID_HANDLE:
        return "DTK error: invalid DTK handle";
    case DTK_UNINITIALIZED:
        return "DTK error: DTK is not initialized";
    case DTK_UNKNOWN:
    default:
        return "DTK error: unknown";
    }
}

} // extern "C"
