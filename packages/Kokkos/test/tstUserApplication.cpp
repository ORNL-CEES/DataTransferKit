/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file   tstUserApplication.cpp
 * \author Stuart Slattery
 * \brief  UserApplication unit tests.
 */
//---------------------------------------------------------------------------//

#include <DTK_ConfigDefs.hpp>
#include <DTK_InputAllocators.hpp>
#include <DTK_ParallelTraits.hpp>
#include <DTK_UserApplication.hpp>
#include <DTK_UserDataInterface.hpp>
#include <DTK_UserFunctionRegistry.hpp>
#include <DTK_View.hpp>

#include <Kokkos_Core.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <memory>

namespace UserAppTest
{
//---------------------------------------------------------------------------//
// User class
//---------------------------------------------------------------------------//
template <class Scalar, class ExecutionSpace>
class UserTestClass
{
  public:
    UserTestClass( const unsigned space_dim, const unsigned size_1,
                   const unsigned size_2, const unsigned offset,
                   const Scalar init_val )
        : _space_dim( space_dim )
        , _size_1( size_1 )
        , _size_2( size_2 )
        , _offset( offset )
        , _data( "test_class_data", size_1, space_dim )
    { /* ... */
    }

  public:
    unsigned _space_dim;
    size_t _size_1;
    size_t _size_2;
    unsigned _offset;
    Kokkos::View<Scalar **> _data;
};

//---------------------------------------------------------------------------//
// User functions.
//---------------------------------------------------------------------------//
// Get the size parameters for building a node list.
template <class Scalar, class ExecutionSpace>
void nodeListSize( std::shared_ptr<void> user_data, unsigned &space_dim,
                   size_t &local_num_nodes, bool &has_ghosts )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_nodes = u->_size_1;
    has_ghosts = true;
}

//---------------------------------------------------------------------------//
// Get the data for a node list.
template <class Scalar, class ExecutionSpace>
void nodeListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> coordinates,
    DataTransferKit::View<bool> is_ghost_node )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t n )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
        {
            coordinates[size_1 * d + n] = n + d + offset;
        }
        is_ghost_node[n] = true;
    };

    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a bounding volume list.
template <class Scalar, class ExecutionSpace>
void boundingVolumeListSize( std::shared_ptr<void> user_data,
                             unsigned &space_dim, size_t &local_num_volumes,
                             bool &has_ghosts )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_volumes = u->_size_1;
    has_ghosts = true;
}

//---------------------------------------------------------------------------//
// Get the data for a bounding volume list.
template <class Scalar, class ExecutionSpace>
void boundingVolumeListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> bounding_volumes,
    DataTransferKit::View<bool> is_ghost_volume )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t v )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
        {
            for ( unsigned h = 0; h < 2; ++h )
            {
                unsigned index = size_1 * space_dim * h + size_1 * d + v;
                bounding_volumes[index] = v + d + h + offset;
            }
        }
        is_ghost_volume[v] = true;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a polyhedron list.
template <class Scalar, class ExecutionSpace>
void polyhedronListSize( std::shared_ptr<void> user_data, unsigned &space_dim,
                         size_t &local_num_nodes, size_t &local_num_faces,
                         size_t &total_nodes_per_face, size_t &local_num_cells,
                         size_t &total_faces_per_cell, bool &has_ghosts )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_nodes = u->_size_1;
    local_num_faces = u->_size_1;
    total_nodes_per_face = u->_size_1;
    local_num_cells = u->_size_1;
    total_faces_per_cell = u->_size_1;
    has_ghosts = true;
}

//---------------------------------------------------------------------------//
// Get the data for a polyhedron list.
template <class Scalar, class ExecutionSpace>
void polyhedronListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> coordinates,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> faces,
    DataTransferKit::View<unsigned> nodes_per_face,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> cells,
    DataTransferKit::View<unsigned> faces_per_cell,
    DataTransferKit::View<int> face_orientation,
    DataTransferKit::View<bool> is_ghost_cell )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t n )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
        {
            coordinates[size_1 * d + n] = n + d + offset;
        }
        faces[n] = n + offset;
        nodes_per_face[n] = n + offset;
        cells[n] = n + offset;
        faces_per_cell[n] = n + offset;
        face_orientation[n] = 1;
        is_ghost_cell[n] = true;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a cell list with a single
// topology.
template <class Scalar, class ExecutionSpace>
void cellListSize( std::shared_ptr<void> user_data, unsigned &space_dim,
                   size_t &local_num_nodes, size_t &local_num_cells,
                   unsigned &nodes_per_cell, bool &has_ghosts )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_nodes = u->_size_1;
    local_num_cells = u->_size_1;
    nodes_per_cell = u->_size_2;
    has_ghosts = true;
}

//---------------------------------------------------------------------------//
// Get the data for a single topology cell list.
template <class Scalar, class ExecutionSpace>
void cellListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> coordinates,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> cells,
    DataTransferKit::View<bool> is_ghost_cell, std::string &cell_topology )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;
    unsigned size_2 = u->_size_2;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t n )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            coordinates[size_1 * d + n] = n + d + offset;
        for ( unsigned v = 0; v < size_2; ++v )
            cells[v * size_1 + n] = n + v + offset;
        is_ghost_cell[n] = true;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();

    cell_topology = "unit_test_topology";
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a cell list with mixed
// topologies.
template <class Scalar, class ExecutionSpace>
void mixedTopologyCellListSize( std::shared_ptr<void> user_data,
                                unsigned &space_dim, size_t &local_num_nodes,
                                size_t &local_num_cells,
                                size_t &total_nodes_per_cell, bool &has_ghosts )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_nodes = u->_size_1;
    local_num_cells = u->_size_1;
    total_nodes_per_cell = u->_size_1;
    has_ghosts = true;
}

//---------------------------------------------------------------------------//
// Get the data for a mixed topology cell list.
template <class Scalar, class ExecutionSpace>
void mixedTopologyCellListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> coordinates,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> cells,
    DataTransferKit::View<unsigned> cell_topology_ids,
    DataTransferKit::View<bool> is_ghost_cell,
    std::vector<std::string> &cell_topologies )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t n )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            coordinates[size_1 * d + n] = n + d + offset;
        cells[n] = n + offset;
        cell_topology_ids[n] = 0;
        is_ghost_cell[n] = true;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();

    cell_topologies.assign( 1, "unit_test_topology" );
}

//---------------------------------------------------------------------------//
// Get the size parameters for a boundary.
template <class Scalar, class ExecutionSpace>
void boundarySize( std::shared_ptr<void> user_data, size_t &local_num_faces )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    local_num_faces = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a boundary.
template <class Scalar, class ExecutionSpace>
void boundaryData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> boundary_cells,
    DataTransferKit::View<unsigned> cell_faces_on_boundary )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned size_1 = u->_size_1;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t n )
    {
        boundary_cells[n] = n + offset;
        cell_faces_on_boundary[n] = n + offset;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Get the size parameters for a degree-of-freedom id map with a single
// number of dofs per object.
template <class Scalar, class ExecutionSpace>
void dofMapSize( std::shared_ptr<void> user_data, size_t &local_num_dofs,
                 size_t &local_num_objects, unsigned &dofs_per_object )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    local_num_dofs = u->_size_1;
    local_num_objects = u->_size_1;
    dofs_per_object = u->_size_2;
}

//---------------------------------------------------------------------------//
// Get the data for a degree-of-freedom id map with a single number of
// dofs per object.
template <class Scalar, class ExecutionSpace>
void dofMapData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::GlobalOrdinal> global_dof_ids,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> object_dof_ids,
    std::string &discretization_type )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned size_1 = u->_size_1;
    unsigned size_2 = u->_size_2;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t n )
    {
        global_dof_ids[n] = n + offset;
        for ( unsigned d = 0; d < size_2; ++d )
            object_dof_ids[size_1 * d + n] = n + d + offset;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();

    discretization_type = "unit_test_discretization";
}

//---------------------------------------------------------------------------//
// Get the size parameters for a degree-of-freedom id map with a
// multiple number of dofs per object (e.g. mixed topology cell lists or
// polyhedron lists).
template <class Scalar, class ExecutionSpace>
void mixedTopologyDOFMapSize( std::shared_ptr<void> user_data,
                              size_t &local_num_dofs, size_t &local_num_objects,
                              size_t &total_dofs_per_object )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    local_num_dofs = u->_size_1;
    local_num_objects = u->_size_1;
    total_dofs_per_object = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a multiple object degree-of-freedom id map
// (e.g. mixed topology cell lists or polyhedron lists).
template <class Scalar, class ExecutionSpace>
void mixedTopologyDOFMapData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::GlobalOrdinal> global_dof_ids,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> object_dof_ids,
    DataTransferKit::View<unsigned> dofs_per_object,
    std::string &discretization_type )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned size_1 = u->_size_1;
    unsigned size_2 = u->_size_2;
    unsigned offset = u->_offset;

    auto fill = KOKKOS_LAMBDA( const size_t n )
    {
        global_dof_ids[n] = n + offset;
        object_dof_ids[n] = n + offset;
        dofs_per_object[n] = size_2;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();

    discretization_type = "unit_test_discretization";
}

//---------------------------------------------------------------------------//
// Get the size parameters for a field. Field must be of size
// local_num_dofs in the associated dof_id_map.
template <class Scalar, class ExecutionSpace>
void fieldSize( std::shared_ptr<void> user_data, unsigned &field_dimension,
                size_t &local_num_dofs )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    field_dimension = u->_space_dim;
    local_num_dofs = u->_size_1;
}

//---------------------------------------------------------------------------//
// Pull data from application into a field.
template <class Scalar, class ExecutionSpace>
void pullFieldData( std::shared_ptr<void> user_data,
                    DataTransferKit::View<Scalar> field_dofs )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;
    auto class_data = u->_data;

    auto pull = KOKKOS_LAMBDA( const size_t n )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            field_dofs[d * size_1 + n] = class_data( n, d );
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          pull );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Push data from a field into the application.
template <class Scalar, class ExecutionSpace>
void pushFieldData( std::shared_ptr<void> user_data,
                    const DataTransferKit::View<Scalar> field_dofs )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;
    auto class_data = u->_data;

    auto push = KOKKOS_LAMBDA( const size_t n )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            class_data( n, d ) = field_dofs[d * size_1 + n];
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          push );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Evaluate a field at a given set of points in a given set of objects.
template <class Scalar, class ExecutionSpace>
void evaluateField(
    std::shared_ptr<void> user_data,
    const DataTransferKit::View<DataTransferKit::Coordinate> evaluation_points,
    const DataTransferKit::View<DataTransferKit::LocalOrdinal> object_ids,
    DataTransferKit::View<Scalar> values )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // The lambda does not properly capture class data so extract it.
    unsigned space_dim = u->_space_dim;
    unsigned size_1 = u->_size_1;

    auto eval = KOKKOS_LAMBDA( const size_t n )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            values[d * size_1 + n] =
                evaluation_points[d * size_1 + n] + object_ids[n];
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          eval );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//

} // end namespace UserAppTest

//---------------------------------------------------------------------------//
// TEST TEMPLATE DECLARATIONS
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, node_list, SC, NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setNodeListSizeFunction(
        UserAppTest::nodeListSize<Scalar, ExecutionSpace>, user_test_class );
    registry->setNodeListDataFunction(
        UserAppTest::nodeListData<Scalar, ExecutionSpace>, user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Get a node list.
    auto node_list = user_app.getNodeList();

    // Check the node list.
    auto host_coordinates = Kokkos::create_mirror_view( node_list.coordinates );
    Kokkos::deep_copy( host_coordinates, node_list.coordinates );
    auto host_is_ghost_node =
        Kokkos::create_mirror_view( node_list.is_ghost_node );
    Kokkos::deep_copy( host_is_ghost_node, node_list.is_ghost_node );
    for ( unsigned i = 0; i < size_1; ++i )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            TEST_EQUALITY( host_coordinates( i, d ), i + d + offset );
        TEST_ASSERT( host_is_ghost_node( i ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, bounding_volume_list, SC,
                                   NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setBoundingVolumeListSizeFunction(
        UserAppTest::boundingVolumeListSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setBoundingVolumeListDataFunction(
        UserAppTest::boundingVolumeListData<Scalar, ExecutionSpace>,
        user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Get a bounding volume list.
    auto bv_list = user_app.getBoundingVolumeList();

    // Check the bounding volumes.
    auto host_bounding_volumes =
        Kokkos::create_mirror_view( bv_list.bounding_volumes );
    Kokkos::deep_copy( host_bounding_volumes, bv_list.bounding_volumes );
    auto host_is_ghost_volume =
        Kokkos::create_mirror_view( bv_list.is_ghost_volume );
    Kokkos::deep_copy( host_is_ghost_volume, bv_list.is_ghost_volume );
    for ( unsigned i = 0; i < size_1; ++i )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            for ( unsigned b = 0; b < 2; ++b )
                TEST_EQUALITY( host_bounding_volumes( i, d, b ),
                               i + d + b + offset );
        TEST_ASSERT( host_is_ghost_volume( i ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, polyhedron_list, SC, NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setPolyhedronListSizeFunction(
        UserAppTest::polyhedronListSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setPolyhedronListDataFunction(
        UserAppTest::polyhedronListData<Scalar, ExecutionSpace>,
        user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Get a polyhedron list.
    auto poly_list = user_app.getPolyhedronList();

    // Check the list.
    auto host_coordinates = Kokkos::create_mirror_view( poly_list.coordinates );
    Kokkos::deep_copy( host_coordinates, poly_list.coordinates );
    auto host_faces = Kokkos::create_mirror_view( poly_list.faces );
    Kokkos::deep_copy( host_faces, poly_list.faces );
    auto host_nodes_per_face =
        Kokkos::create_mirror_view( poly_list.nodes_per_face );
    Kokkos::deep_copy( host_nodes_per_face, poly_list.nodes_per_face );
    auto host_cells = Kokkos::create_mirror_view( poly_list.cells );
    Kokkos::deep_copy( host_cells, poly_list.cells );
    auto host_faces_per_cell =
        Kokkos::create_mirror_view( poly_list.faces_per_cell );
    Kokkos::deep_copy( host_faces_per_cell, poly_list.faces_per_cell );
    auto host_face_orientation =
        Kokkos::create_mirror_view( poly_list.face_orientation );
    Kokkos::deep_copy( host_face_orientation, poly_list.face_orientation );
    auto host_is_ghost_cell =
        Kokkos::create_mirror_view( poly_list.is_ghost_cell );
    Kokkos::deep_copy( host_is_ghost_cell, poly_list.is_ghost_cell );
    for ( unsigned i = 0; i < size_1; ++i )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            TEST_EQUALITY( host_coordinates( i, d ), i + d + offset );
        TEST_EQUALITY( host_faces( i ), i + offset );
        TEST_EQUALITY( host_nodes_per_face( i ), i + offset );
        TEST_EQUALITY( host_cells( i ), i + offset );
        TEST_EQUALITY( host_faces_per_cell( i ), i + offset );
        TEST_EQUALITY( host_face_orientation( i ), 1 );
        TEST_ASSERT( host_is_ghost_cell( i ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, single_topology_cell, SC,
                                   NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setCellListSizeFunction(
        UserAppTest::cellListSize<Scalar, ExecutionSpace>, user_test_class );
    registry->setCellListDataFunction(
        UserAppTest::cellListData<Scalar, ExecutionSpace>, user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Get a cell list.
    std::vector<std::string> cell_topologies;
    auto cell_list = user_app.getCellList( cell_topologies );
    TEST_EQUALITY( cell_list.cells.rank(), 2 );

    // Check the list.
    auto host_coordinates = Kokkos::create_mirror_view( cell_list.coordinates );
    Kokkos::deep_copy( host_coordinates, cell_list.coordinates );
    auto host_cells = Kokkos::create_mirror_view( cell_list.cells );
    Kokkos::deep_copy( host_cells, cell_list.cells );
    auto host_is_ghost_cell =
        Kokkos::create_mirror_view( cell_list.is_ghost_cell );
    Kokkos::deep_copy( host_is_ghost_cell, cell_list.is_ghost_cell );
    for ( unsigned i = 0; i < size_1; ++i )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            TEST_EQUALITY( host_coordinates( i, d ), i + d + offset );
        for ( unsigned v = 0; v < size_2; ++v )
            TEST_EQUALITY( host_cells( i, v ), i + v + offset );
        TEST_ASSERT( host_is_ghost_cell( i ) );
    }
    TEST_EQUALITY( cell_topologies.size(), 1 );
    TEST_EQUALITY( cell_topologies[0], "unit_test_topology" );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, multiple_topology_cell, SC,
                                   NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setMixedTopologyCellListSizeFunction(
        UserAppTest::mixedTopologyCellListSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setMixedTopologyCellListDataFunction(
        UserAppTest::mixedTopologyCellListData<Scalar, ExecutionSpace>,
        user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Get a cell list.
    std::vector<std::string> cell_topologies;
    auto cell_list = user_app.getCellList( cell_topologies );
    TEST_EQUALITY( cell_list.cells.rank(), 1 );

    // Check the list.
    auto host_coordinates = Kokkos::create_mirror_view( cell_list.coordinates );
    Kokkos::deep_copy( host_coordinates, cell_list.coordinates );
    auto host_cells = Kokkos::create_mirror_view( cell_list.cells );
    Kokkos::deep_copy( host_cells, cell_list.cells );
    auto host_cell_topology_ids =
        Kokkos::create_mirror_view( cell_list.cell_topology_ids );
    Kokkos::deep_copy( host_cell_topology_ids, cell_list.cell_topology_ids );
    auto host_is_ghost_cell =
        Kokkos::create_mirror_view( cell_list.is_ghost_cell );
    Kokkos::deep_copy( host_is_ghost_cell, cell_list.is_ghost_cell );
    for ( unsigned i = 0; i < size_1; ++i )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            TEST_EQUALITY( host_coordinates( i, d ), i + d + offset );
        TEST_EQUALITY( host_cells( i ), i + offset );
        TEST_EQUALITY( host_cell_topology_ids( i ), 0 );
        TEST_ASSERT( host_is_ghost_cell( i ) );
    }
    TEST_EQUALITY( cell_topologies.size(), 1 );
    TEST_EQUALITY( cell_topologies[0], "unit_test_topology" );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, boundary, SC, NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    std::string boundary_name( "unit_test_boundary" );
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setBoundarySizeFunction(
        boundary_name, UserAppTest::boundarySize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setBoundaryDataFunction(
        boundary_name, UserAppTest::boundaryData<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setCellListSizeFunction(
        UserAppTest::cellListSize<Scalar, ExecutionSpace>, user_test_class );
    registry->setCellListDataFunction(
        UserAppTest::cellListData<Scalar, ExecutionSpace>, user_test_class );
    registry->setPolyhedronListSizeFunction(
        UserAppTest::polyhedronListSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setPolyhedronListDataFunction(
        UserAppTest::polyhedronListData<Scalar, ExecutionSpace>,
        user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Test with a cell list.
    {
        // Create a cell list.
        std::vector<std::string> discretization;
        auto cell_list = user_app.getCellList( discretization );

        // Get the boundary of the list.
        user_app.getBoundary( boundary_name, cell_list );

        // Check the boundary.
        auto host_boundary_cells =
            Kokkos::create_mirror_view( cell_list.boundary_cells );
        Kokkos::deep_copy( host_boundary_cells, cell_list.boundary_cells );
        auto host_cell_faces_on_boundary =
            Kokkos::create_mirror_view( cell_list.cell_faces_on_boundary );
        Kokkos::deep_copy( host_cell_faces_on_boundary,
                           cell_list.cell_faces_on_boundary );
        for ( unsigned i = 0; i < size_1; ++i )
        {
            TEST_EQUALITY( host_boundary_cells( i ), i + offset );
            TEST_EQUALITY( host_cell_faces_on_boundary( i ), i + offset );
        }
    }

    // Test with a cell list.
    {
        // Create a polyhedron list.
        auto poly_list = user_app.getPolyhedronList();

        // Get the boundary of the list.
        user_app.getBoundary( boundary_name, poly_list );

        // Check the boundary.
        auto host_boundary_cells =
            Kokkos::create_mirror_view( poly_list.boundary_cells );
        Kokkos::deep_copy( host_boundary_cells, poly_list.boundary_cells );
        auto host_cell_faces_on_boundary =
            Kokkos::create_mirror_view( poly_list.cell_faces_on_boundary );
        Kokkos::deep_copy( host_cell_faces_on_boundary,
                           poly_list.cell_faces_on_boundary );
        for ( unsigned i = 0; i < size_1; ++i )
        {
            TEST_EQUALITY( host_boundary_cells( i ), i + offset );
            TEST_EQUALITY( host_cell_faces_on_boundary( i ), i + offset );
        }
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, single_topology_dof, SC,
                                   NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setDOFMapSizeFunction(
        UserAppTest::dofMapSize<Scalar, ExecutionSpace>, user_test_class );
    registry->setDOFMapDataFunction(
        UserAppTest::dofMapData<Scalar, ExecutionSpace>, user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Create a map.
    std::string discretization_type;
    auto dof_map = user_app.getDOFMap( discretization_type );

    // Check the map.
    TEST_EQUALITY( dof_map.object_dof_ids.rank(), 2 );
    auto host_global_dof_ids =
        Kokkos::create_mirror_view( dof_map.global_dof_ids );
    Kokkos::deep_copy( host_global_dof_ids, dof_map.global_dof_ids );
    auto host_object_dof_ids =
        Kokkos::create_mirror_view( dof_map.object_dof_ids );
    Kokkos::deep_copy( host_object_dof_ids, dof_map.object_dof_ids );
    for ( unsigned i = 0; i < size_1; ++i )
    {
        host_global_dof_ids( i ) = i + offset;
        for ( unsigned d = 0; d < size_2; ++d )
            host_object_dof_ids( i, d ) = i + d + offset;
    }
    TEST_EQUALITY( discretization_type, "unit_test_discretization" );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, multiple_topology_dof, SC,
                                   NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );
    registry->setMixedTopologyDOFMapSizeFunction(
        UserAppTest::mixedTopologyDOFMapSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setMixedTopologyDOFMapDataFunction(
        UserAppTest::mixedTopologyDOFMapData<Scalar, ExecutionSpace>,
        user_test_class );

    // Create a map.
    std::string discretization_type;
    auto dof_map = user_app.getDOFMap( discretization_type );

    // Check the map.
    TEST_EQUALITY( dof_map.object_dof_ids.rank(), 1 );
    auto host_global_dof_ids =
        Kokkos::create_mirror_view( dof_map.global_dof_ids );
    Kokkos::deep_copy( host_global_dof_ids, dof_map.global_dof_ids );
    auto host_object_dof_ids =
        Kokkos::create_mirror_view( dof_map.object_dof_ids );
    Kokkos::deep_copy( host_object_dof_ids, dof_map.object_dof_ids );
    auto host_dofs_per_object =
        Kokkos::create_mirror_view( dof_map.dofs_per_object );
    Kokkos::deep_copy( host_dofs_per_object, dof_map.dofs_per_object );
    for ( unsigned i = 0; i < size_1; ++i )
    {
        host_global_dof_ids( i ) = i + offset;
        host_object_dof_ids( i ) = i + offset;
        host_dofs_per_object( i ) = size_2;
    }
    TEST_EQUALITY( discretization_type, "unit_test_discretization" );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, field_push_pull, SC, NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    std::string field_name = "test_field";
    registry->setFieldSizeFunction(
        field_name, UserAppTest::fieldSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setPullFieldDataFunction(
        field_name, UserAppTest::pullFieldData<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setPushFieldDataFunction(
        field_name, UserAppTest::pushFieldData<Scalar, ExecutionSpace>,
        user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Create a field.
    auto field_1 = user_app.getField( field_name );

    // Put some data in the field.
    auto fill_field = KOKKOS_LAMBDA( const size_t i )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            field_1.dofs( i, d ) = i + d;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill_field );
    Kokkos::fence();

    // Push the field into the app.
    user_app.pushField( field_name, field_1 );

    // Create a second field.
    auto field_2 = user_app.getField( field_name );

    // Pull the field out of the app.
    user_app.pullField( field_name, field_2 );

    // Check the pulled field.
    auto host_dofs = Kokkos::create_mirror_view( field_2.dofs );
    Kokkos::deep_copy( host_dofs, field_2.dofs );
    for ( unsigned i = 0; i < size_1; ++i )
        for ( unsigned d = 0; d < space_dim; ++d )
            TEST_EQUALITY( host_dofs( i, d ), i + d );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, field_eval, SC, NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    std::string field_name = "test_field";
    registry->setFieldSizeFunction(
        field_name, UserAppTest::fieldSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setEvaluateFieldFunction(
        field_name, UserAppTest::evaluateField<Scalar, ExecutionSpace>,
        user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Create an evaluation set.
    auto eval_set = DataTransferKit::InputAllocators<
        Kokkos::LayoutLeft, ExecutionSpace>::allocateEvaluationSet( size_1,
                                                                    space_dim );
    auto fill_eval_set = KOKKOS_LAMBDA( const size_t i )
    {
        for ( unsigned d = 0; d < space_dim; ++d )
            eval_set.evaluation_points( i, d ) = i + d;
        eval_set.object_ids( i ) = i;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill_eval_set );
    Kokkos::fence();

    // Create a field.
    auto field = user_app.getField( field_name );

    // Evaluate the field.
    user_app.evaluateField( field_name, eval_set, field );

    // Check the evaluation.
    auto host_dofs = Kokkos::create_mirror_view( field.dofs );
    Kokkos::deep_copy( host_dofs, field.dofs );
    auto host_points = Kokkos::create_mirror_view( eval_set.evaluation_points );
    Kokkos::deep_copy( host_points, eval_set.evaluation_points );
    auto host_object_ids = Kokkos::create_mirror_view( eval_set.object_ids );
    Kokkos::deep_copy( host_object_ids, eval_set.object_ids );
    for ( unsigned i = 0; i < size_1; ++i )
        for ( unsigned d = 0; d < space_dim; ++d )
            TEST_EQUALITY( host_dofs( i, d ),
                           host_points( i, d ) + host_object_ids( i ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, missing_function, SC, NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions. Forget to set the data function on
    // purpose.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setNodeListSizeFunction(
        UserAppTest::nodeListSize<Scalar, ExecutionSpace>, user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // Get a node list. This should throw because the function is missing.
    bool caught_exception = false;
    try
    {
        auto node_list = user_app.getNodeList();
    }
    catch ( DataTransferKit::DataTransferKitException &e )
    {
        caught_exception = true;
    }
    TEST_ASSERT( caught_exception );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, too_many_functions, SC, NO )
{
    // Test types.
    using DeviceType = typename NO::device_type;
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    const unsigned space_dim = 3;
    const unsigned size_1 = 100;
    const unsigned size_2 = 5;
    const unsigned offset = 8;
    const SC init_val = Teuchos::ScalarTraits<SC>::one();
    auto user_test_class =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>(
            space_dim, size_1, size_2, offset, init_val );

    // Set the user functions. Set both single and mixed topology
    // functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setCellListSizeFunction(
        UserAppTest::cellListSize<Scalar, ExecutionSpace>, user_test_class );
    registry->setCellListDataFunction(
        UserAppTest::cellListData<Scalar, ExecutionSpace>, user_test_class );
    registry->setDOFMapSizeFunction(
        UserAppTest::dofMapSize<Scalar, ExecutionSpace>, user_test_class );
    registry->setDOFMapDataFunction(
        UserAppTest::dofMapData<Scalar, ExecutionSpace>, user_test_class );
    registry->setMixedTopologyCellListSizeFunction(
        UserAppTest::mixedTopologyCellListSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setMixedTopologyCellListDataFunction(
        UserAppTest::mixedTopologyCellListData<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setMixedTopologyDOFMapSizeFunction(
        UserAppTest::mixedTopologyDOFMapSize<Scalar, ExecutionSpace>,
        user_test_class );
    registry->setMixedTopologyDOFMapDataFunction(
        UserAppTest::mixedTopologyDOFMapData<Scalar, ExecutionSpace>,
        user_test_class );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    // First get a cell list. We registered both mixed and single topology
    // function so this will fail.
    bool caught_exception = false;
    try
    {
        std::vector<std::string> cell_topologies;
        auto cell_list = user_app.getCellList( cell_topologies );
    }
    catch ( DataTransferKit::DataTransferKitException &e )
    {
        caught_exception = true;
    }
    TEST_ASSERT( caught_exception );

    // Next get a dof id map. We registered both mixed and single topology
    // function so this will fail.
    caught_exception = false;
    try
    {
        std::string discretization_type;
        auto dof_map = user_app.getDOFMap( discretization_type );
    }
    catch ( DataTransferKit::DataTransferKitException &e )
    {
        caught_exception = true;
    }
    TEST_ASSERT( caught_exception );
}

//---------------------------------------------------------------------------//
// TEST TEMPLATE INSTANTIATIONS
//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( SCALAR, NODE )                                        \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, node_list, SCALAR,  \
                                          NODE )                               \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication,                     \
                                          bounding_volume_list, SCALAR, NODE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, polyhedron_list,    \
                                          SCALAR, NODE )                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication,                     \
                                          single_topology_cell, SCALAR, NODE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                      \
        UserApplication, multiple_topology_cell, SCALAR, NODE )                \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, boundary, SCALAR,   \
                                          NODE )                               \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication,                     \
                                          single_topology_dof, SCALAR, NODE )  \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                      \
        UserApplication, multiple_topology_dof, SCALAR, NODE )                 \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, field_push_pull,    \
                                          SCALAR, NODE )                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, field_eval, SCALAR, \
                                          NODE )                               \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, missing_function,   \
                                          SCALAR, NODE )                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, too_many_functions, \
                                          SCALAR, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SN( UNIT_TEST_GROUP )

//---------------------------------------------------------------------------//
// end tstUserApplication.cpp
//---------------------------------------------------------------------------//
