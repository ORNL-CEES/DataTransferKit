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

#include <DTK_CellTypes.h>
#include <DTK_ConfigDefs.hpp>
#include <DTK_InputAllocators.hpp>
#include <DTK_ParallelTraits.hpp>
#include <DTK_UserApplication.hpp>
#include <DTK_UserDataInterface.hpp>
#include <DTK_UserFunctionRegistry.hpp>
#include <DTK_View.hpp>

#include "DTK_TestApplicationHelpers.hpp"

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
struct UserTestClass
{
    UserTestClass()
        : _data( "test_class_data", _size_1, _space_dim )
    { /* ... */
    }

    const unsigned _space_dim = 3;
    const size_t _size_1 = 100;
    const size_t _size_2 = 5;
    const unsigned _offset = 8;
    const std::string _boundary_name = "unit_test_boundary";
    const std::string _field_name = "test_field";
    Kokkos::View<Scalar **> _data;
};

//---------------------------------------------------------------------------//
// User functions.
//---------------------------------------------------------------------------//
// Get the size parameters for building a node list.
template <class Scalar, class ExecutionSpace>
void nodeListSize( std::shared_ptr<void> user_data, unsigned &space_dim,
                   size_t &local_num_nodes )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_nodes = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a node list.
template <class Scalar, class ExecutionSpace>
void nodeListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> coordinates )
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
    };

    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a bounding volume list.
template <class Scalar, class ExecutionSpace>
void boundingVolumeListSize( std::shared_ptr<void> user_data,
                             unsigned &space_dim, size_t &local_num_volumes )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_volumes = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a bounding volume list.
template <class Scalar, class ExecutionSpace>
void boundingVolumeListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> bounding_volumes )
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
                         size_t &total_face_nodes, size_t &local_num_cells,
                         size_t &total_cell_faces )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_nodes = u->_size_1;
    local_num_faces = u->_size_1;
    total_face_nodes = u->_size_1;
    local_num_cells = u->_size_1;
    total_cell_faces = u->_size_1;
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
    DataTransferKit::View<int> face_orientation )
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
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Get the size parameters for building a cell list.
template <class Scalar, class ExecutionSpace>
void cellListSize( std::shared_ptr<void> user_data, unsigned &space_dim,
                   size_t &local_num_nodes, size_t &local_num_cells,
                   size_t &total_cell_nodes )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    space_dim = u->_space_dim;
    local_num_nodes = u->_size_1;
    local_num_cells = u->_size_1;
    total_cell_nodes = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a cell list.
template <class Scalar, class ExecutionSpace>
void cellListData(
    std::shared_ptr<void> user_data,
    DataTransferKit::View<DataTransferKit::Coordinate> coordinates,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> cells,
    DataTransferKit::View<DTK_CellTopology> cell_topologies )
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
        cell_topologies[n] = DTK_TET_4;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, size_1 ),
                          fill );
    Kokkos::fence();
}

//---------------------------------------------------------------------------//
// Get the size parameters for a boundary.
template <class Scalar, class ExecutionSpace>
void boundarySize( std::shared_ptr<void> user_data,
                   const std::string &boundary_name, size_t &local_num_faces )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // Here one could do actions depening on the name, but in the tests we
    // simply ignore it
    (void)boundary_name;
    local_num_faces = u->_size_1;
}

//---------------------------------------------------------------------------//
// Get the data for a boundary.
template <class Scalar, class ExecutionSpace>
void boundaryData(
    std::shared_ptr<void> user_data, const std::string &boundary_name,
    DataTransferKit::View<DataTransferKit::LocalOrdinal> boundary_cells,
    DataTransferKit::View<unsigned> cell_faces_on_boundary )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // Here one could do actions depening on the name, but in the tests we
    // simply ignore it
    (void)boundary_name;

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
void fieldSize( std::shared_ptr<void> user_data, const std::string &field_name,
                unsigned &field_dimension, size_t &local_num_dofs )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // Here one could do actions depening on the name, but in the tests we
    // simply ignore it
    (void)field_name;

    field_dimension = u->_space_dim;
    local_num_dofs = u->_size_1;
}

//---------------------------------------------------------------------------//
// Pull data from application into a field.
template <class Scalar, class ExecutionSpace>
void pullFieldData( std::shared_ptr<void> user_data,
                    const std::string &field_name,
                    DataTransferKit::View<Scalar> field_dofs )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // Here one could do actions depening on the name, but in the tests we
    // simply ignore it
    (void)field_name;

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
                    const std::string &field_name,
                    const DataTransferKit::View<Scalar> field_dofs )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // Here one could do actions depening on the name, but in the tests we
    // simply ignore it
    (void)field_name;

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
    std::shared_ptr<void> user_data, const std::string &field_name,
    const DataTransferKit::View<DataTransferKit::Coordinate> evaluation_points,
    const DataTransferKit::View<DataTransferKit::LocalOrdinal> object_ids,
    DataTransferKit::View<Scalar> values )
{
    auto u = std::static_pointer_cast<UserTestClass<Scalar, ExecutionSpace>>(
        user_data );

    // Here one could do actions depening on the name, but in the tests we
    // simply ignore it
    (void)field_name;

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
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, node_list, SC, DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setNodeListSizeFunction(
        UserAppTest::nodeListSize<Scalar, ExecutionSpace>, u );
    registry->setNodeListDataFunction(
        UserAppTest::nodeListData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_node_list( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, bounding_volume_list, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setBoundingVolumeListSizeFunction(
        UserAppTest::boundingVolumeListSize<Scalar, ExecutionSpace>, u );
    registry->setBoundingVolumeListDataFunction(
        UserAppTest::boundingVolumeListData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_bounding_volume_list( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, polyhedron_list, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setPolyhedronListSizeFunction(
        UserAppTest::polyhedronListSize<Scalar, ExecutionSpace>, u );
    registry->setPolyhedronListDataFunction(
        UserAppTest::polyhedronListData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_polyhedron_list( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, multiple_topology_cell, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setCellListSizeFunction(
        UserAppTest::cellListSize<Scalar, ExecutionSpace>, u );
    registry->setCellListDataFunction(
        UserAppTest::cellListData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_multiple_topology_cell( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, boundary, SC, DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setBoundarySizeFunction(
        UserAppTest::boundarySize<Scalar, ExecutionSpace>, u );
    registry->setBoundaryDataFunction(
        UserAppTest::boundaryData<Scalar, ExecutionSpace>, u );
    registry->setCellListSizeFunction(
        UserAppTest::cellListSize<Scalar, ExecutionSpace>, u );
    registry->setCellListDataFunction(
        UserAppTest::cellListData<Scalar, ExecutionSpace>, u );
    registry->setPolyhedronListSizeFunction(
        UserAppTest::polyhedronListSize<Scalar, ExecutionSpace>, u );
    registry->setPolyhedronListDataFunction(
        UserAppTest::polyhedronListData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_boundary( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, single_topology_dof, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setDOFMapSizeFunction(
        UserAppTest::dofMapSize<Scalar, ExecutionSpace>, u );
    registry->setDOFMapDataFunction(
        UserAppTest::dofMapData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_single_topology_dof( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, multiple_topology_dof, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setMixedTopologyDOFMapSizeFunction(
        UserAppTest::mixedTopologyDOFMapSize<Scalar, ExecutionSpace>, u );
    registry->setMixedTopologyDOFMapDataFunction(
        UserAppTest::mixedTopologyDOFMapData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_multiple_topology_dof( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, field_push_pull, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setFieldSizeFunction(
        UserAppTest::fieldSize<Scalar, ExecutionSpace>, u );
    registry->setPullFieldDataFunction(
        UserAppTest::pullFieldData<Scalar, ExecutionSpace>, u );
    registry->setPushFieldDataFunction(
        UserAppTest::pushFieldData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_field_push_pull( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, field_eval, SC, DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setFieldSizeFunction(
        UserAppTest::fieldSize<Scalar, ExecutionSpace>, u );
    registry->setEvaluateFieldFunction(
        UserAppTest::evaluateField<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_field_eval( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, missing_function, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions. Forget to set the data function on
    // purpose.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setNodeListSizeFunction(
        UserAppTest::nodeListSize<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_missing_function( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( UserApplication, too_many_functions, SC,
                                   DeviceType )
{
    // Test types.
    using ExecutionSpace = typename DeviceType::execution_space;
    using Scalar = SC;

    // Create the test class.
    auto u =
        std::make_shared<UserAppTest::UserTestClass<Scalar, ExecutionSpace>>();

    // Set the user functions. Set both single and mixed topology
    // functions.
    auto registry =
        std::make_shared<DataTransferKit::UserFunctionRegistry<Scalar>>();
    registry->setDOFMapSizeFunction(
        UserAppTest::dofMapSize<Scalar, ExecutionSpace>, u );
    registry->setDOFMapDataFunction(
        UserAppTest::dofMapData<Scalar, ExecutionSpace>, u );
    registry->setMixedTopologyDOFMapSizeFunction(
        UserAppTest::mixedTopologyDOFMapSize<Scalar, ExecutionSpace>, u );
    registry->setMixedTopologyDOFMapDataFunction(
        UserAppTest::mixedTopologyDOFMapData<Scalar, ExecutionSpace>, u );

    // Create the user application.
    DataTransferKit::UserApplication<Scalar, ExecutionSpace> user_app(
        registry );

    test_too_many_functions( user_app, *u, out, success );
}

//---------------------------------------------------------------------------//
// TEST TEMPLATE INSTANTIATIONS
//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitInterface_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( SCALAR, NODE )                                        \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, node_list, SCALAR,  \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                      \
        UserApplication, bounding_volume_list, SCALAR, DeviceType##NODE )      \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, polyhedron_list,    \
                                          SCALAR, DeviceType##NODE )           \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                      \
        UserApplication, multiple_topology_cell, SCALAR, DeviceType##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, boundary, SCALAR,   \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                      \
        UserApplication, single_topology_dof, SCALAR, DeviceType##NODE )       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                      \
        UserApplication, multiple_topology_dof, SCALAR, DeviceType##NODE )     \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, field_push_pull,    \
                                          SCALAR, DeviceType##NODE )           \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, field_eval, SCALAR, \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, missing_function,   \
                                          SCALAR, DeviceType##NODE )           \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( UserApplication, too_many_functions, \
                                          SCALAR, DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SN( UNIT_TEST_GROUP )

//---------------------------------------------------------------------------//
// end tstUserApplication.cpp
//---------------------------------------------------------------------------//
