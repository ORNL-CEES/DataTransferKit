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
 * \file   tstMapInterface.cpp
 * \author Stuart Slattery
 * \brief  DTK view unit tests.
 */
//---------------------------------------------------------------------------//

#include <DTK_C_API.h>
#include <DTK_ParallelTraits.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Kokkos_Core.hpp>

#include <memory>

//---------------------------------------------------------------------------//
// User implementation
template <class Space>
struct TestUserData
{
    Kokkos::View<double * [3], Space> coords;
    Kokkos::View<double *, Space> field;

    TestUserData( const int size )
        : coords( "coords", size )
        , field( "field", size )
    {
    }
};

template <class Space>
void nodeListSize( void *user_data, unsigned *space_dim,
                   size_t *local_num_nodes )
{
    TestUserData<Space> *data = static_cast<TestUserData<Space> *>( user_data );
    *space_dim = data->coords.extent( 1 );
    *local_num_nodes = data->coords.extent( 0 );
}

template <class Space>
void nodeListData( void *user_data, Coordinate *coords )
{
    TestUserData<Space> *data = static_cast<TestUserData<Space> *>( user_data );
    int num_node = data->coords.extent( 0 );
    for ( unsigned n = 0; n < data->coords.extent( 0 ); ++n )
        for ( unsigned d = 0; d < data->coords.extent( 1 ); ++d )
            coords[num_node * d + n] = data->coords( n, d );
}

template <class Space>
void fieldSize( void *user_data, const char *field_name,
                unsigned *field_dimension, size_t *local_num_dofs )
{
    TestUserData<Space> *data = static_cast<TestUserData<Space> *>( user_data );
    *field_dimension = 1;
    *local_num_dofs = data->field.extent( 0 );
}

template <class Space>
void pullField( void *user_data, const char *, double *field_dofs )
{
    TestUserData<Space> *data = static_cast<TestUserData<Space> *>( user_data );
    for ( unsigned i = 0; i < data->field.extent( 0 ); ++i )
        field_dofs[i] = data->field( i );
}
template <class Space>
void pushField( void *user_data, const char *, const double *field_dofs )
{
    TestUserData<Space> *data = static_cast<TestUserData<Space> *>( user_data );
    for ( unsigned i = 0; i < data->field.extent( 0 ); ++i )
        data->field( i ) = field_dofs[i];
}

//---------------------------------------------------------------------------//
// Test execution space enumeration selector.
template <class Space>
struct SpaceSelector;

#if defined( KOKKOS_ENABLE_SERIAL )
template <>
struct SpaceSelector<DataTransferKit::Serial>
{
    static constexpr DTK_ExecutionSpace value() { return DTK_SERIAL; }
};
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
template <>
struct SpaceSelector<DataTransferKit::OpenMP>
{
    static constexpr DTK_ExecutionSpace value() { return DTK_OPENMP; }
};
#endif

#if defined( KOKKOS_ENABLE_SERIAL ) || defined( KOKKOS_ENABLE_OPENMP )
template <>
struct SpaceSelector<DataTransferKit::HostSpace>
{
    static constexpr DTK_MemorySpace value() { return DTK_HOST_SPACE; }
};
#endif

#if defined( KOKKOS_ENABLE_CUDA )
template <>
struct SpaceSelector<DataTransferKit::Cuda>
{
    static constexpr DTK_ExecutionSpace value() { return DTK_CUDA; }
};

template <>
struct SpaceSelector<DataTransferKit::CudaUVMSpace>
{
    static constexpr DTK_MemorySpace value() { return DTK_CUDAUVM_SPACE; }
};
#endif

//---------------------------------------------------------------------------//
// Run the test.
template <class MapSpace, class SourceSpace, class TargetSpace>
void test( bool &success, Teuchos::FancyOStream &out )
{
    // Initialize DTK. The test harness initializes kokkos already.
    DTK_initialize();
    TEST_EQUALITY( errno, DTK_SUCCESS );

    // Check error handling on a bad map.
    DTK_MapHandle bad_handle = nullptr;
    DTK_applyMap( bad_handle, "bad", "bad" );
    TEST_EQUALITY( errno, DTK_INVALID_HANDLE );
    DTK_destroyMap( bad_handle );
    TEST_EQUALITY( errno, DTK_INVALID_HANDLE );

    // Get the communicator
    auto teuchos_comm = Teuchos::DefaultComm<int>::getComm();

    // Create a source and target data set.
    int num_point = 1000;
    auto src_data = std::make_shared<TestUserData<SourceSpace>>( num_point );
    auto tgt_data = std::make_shared<TestUserData<TargetSpace>>( num_point );

    // Assign source and target points and some field data. Make the points
    // the same for now to test the nearest neighbor operator. Invert the rank
    // though so at least they live on different processors.
    int comm_rank = teuchos_comm->getRank();
    int inverse_rank = teuchos_comm->getSize() - comm_rank - 1;
    for ( int p = 0; p < num_point; ++p )
    {
        for ( int d = 0; d < 3; ++d )
        {

            src_data->coords( p, d ) = 1.0 * p + comm_rank * num_point;
            tgt_data->coords( p, d ) = 1.0 * p + inverse_rank * num_point;
        }
        src_data->field( p ) = 1.0 * p + comm_rank * num_point;
        tgt_data->field( p ) = 0.0;
    }

    // Create the source user application instance.
    auto src_handle =
        DTK_createUserApplication( SpaceSelector<SourceSpace>::value() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( src_handle, DTK_NODE_LIST_SIZE_FUNCTION,
                         ( void ( * )() ) & nodeListSize<SourceSpace>,
                         src_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( src_handle, DTK_NODE_LIST_DATA_FUNCTION,
                         ( void ( * )() ) & nodeListData<SourceSpace>,
                         src_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( src_handle, DTK_FIELD_SIZE_FUNCTION,
                         ( void ( * )() ) & fieldSize<SourceSpace>,
                         src_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( src_handle, DTK_PULL_FIELD_DATA_FUNCTION,
                         ( void ( * )() ) & pullField<SourceSpace>,
                         src_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );

    // Create the target user application instance.
    auto tgt_handle =
        DTK_createUserApplication( SpaceSelector<TargetSpace>::value() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( tgt_handle, DTK_NODE_LIST_SIZE_FUNCTION,
                         ( void ( * )() ) & nodeListSize<TargetSpace>,
                         tgt_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( tgt_handle, DTK_NODE_LIST_DATA_FUNCTION,
                         ( void ( * )() ) & nodeListData<TargetSpace>,
                         tgt_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( tgt_handle, DTK_FIELD_SIZE_FUNCTION,
                         ( void ( * )() ) & fieldSize<TargetSpace>,
                         tgt_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_setUserFunction( tgt_handle, DTK_PUSH_FIELD_DATA_FUNCTION,
                         ( void ( * )() ) & pushField<TargetSpace>,
                         tgt_data.get() );
    TEST_EQUALITY( errno, DTK_SUCCESS );

    // Create a map.
    std::string const options = R"({ "Map Type": "Nearest Neighbor" })";
    auto comm = Teuchos::getRawMpiComm( *teuchos_comm );
    auto map_handle = DTK_createMap( SpaceSelector<MapSpace>::value(), comm,
                                     src_handle, tgt_handle, options.c_str() );
    TEST_EQUALITY( errno, DTK_SUCCESS );

    // Apply the map.
    DTK_applyMap( map_handle, "dummy", "dummy" );
    TEST_EQUALITY( errno, DTK_SUCCESS );

    // Check the results.
    double const relative_tolerance = 1e-14;
    // NOTE adding the same value to both lhs and rhs to resolve floating point
    // comparison issues with zero using Teuchos assertion macro
    double const shift_from_zero = 3.14;
    for ( int p = 0; p < num_point; ++p )
    {
        TEST_FLOATING_EQUALITY( tgt_data->field( p ) + shift_from_zero,
                                1.0 * p + inverse_rank * num_point +
                                    shift_from_zero,
                                relative_tolerance );
    }

    // Cleanup.
    DTK_destroyMap( map_handle );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_destroyUserApplication( src_handle );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_destroyUserApplication( tgt_handle );
    TEST_EQUALITY( errno, DTK_SUCCESS );
    DTK_finalize();
    TEST_EQUALITY( errno, DTK_SUCCESS );
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
#if defined( KOKKOS_ENABLE_SERIAL )
TEUCHOS_UNIT_TEST( MapInterface, Serial )
{
    test<DataTransferKit::Serial, DataTransferKit::HostSpace,
         DataTransferKit::HostSpace>( success, out );
}
#endif

//---------------------------------------------------------------------------//
#if defined( KOKKOS_ENABLE_OPENMP )
TEUCHOS_UNIT_TEST( MapInterface, OpenMP )
{
    test<DataTransferKit::OpenMP, DataTransferKit::HostSpace,
         DataTransferKit::HostSpace>( success, out );
}
#endif

//---------------------------------------------------------------------------//
#if defined( KOKKOS_ENABLE_CUDA )
TEUCHOS_UNIT_TEST( MapInterface, Cuda )
{
    test<DataTransferKit::Cuda, DataTransferKit::CudaUVMSpace,
         DataTransferKit::CudaUVMSpace>( success, out );
}
#endif

//---------------------------------------------------------------------------//
#if defined( KOKKOS_ENABLE_SERIAL ) && defined( KOKKOS_ENABLE_CUDA )
TEUCHOS_UNIT_TEST( MapInterface, SerialAndCuda )
{
    test<DataTransferKit::Serial, DataTransferKit::CudaUVMSpace,
         DataTransferKit::CudaUVMSpace>( success, out );

    test<DataTransferKit::Serial, DataTransferKit::HostSpace,
         DataTransferKit::CudaUVMSpace>( success, out );

    test<DataTransferKit::Serial, DataTransferKit::CudaUVMSpace,
         DataTransferKit::HostSpace>( success, out );
}
#endif

//---------------------------------------------------------------------------//
#if defined( KOKKOS_ENABLE_OPENMP ) && defined( KOKKOS_ENABLE_CUDA )
TEUCHOS_UNIT_TEST( MapInterface, OpenmpAndCuda )
{
    test<DataTransferKit::OpenMP, DataTransferKit::CudaUVMSpace,
         DataTransferKit::CudaUVMSpace>( success, out );

    test<DataTransferKit::OpenMP, DataTransferKit::HostSpace,
         DataTransferKit::CudaUVMSpace>( success, out );

    test<DataTransferKit::OpenMP, DataTransferKit::CudaUVMSpace,
         DataTransferKit::HostSpace>( success, out );
}
#endif

//---------------------------------------------------------------------------//
#if ( defined( KOKKOS_ENABLE_SERIAL ) || defined( KOKKOS_ENABLE_OPENMP ) ) &&  \
    defined( KOKKOS_ENABLE_CUDA )
TEUCHOS_UNIT_TEST( MapInterface, HostAndCuda )
{
    test<DataTransferKit::Cuda, DataTransferKit::HostSpace,
         DataTransferKit::CudaUVMSpace>( success, out );

    test<DataTransferKit::Cuda, DataTransferKit::CudaUVMSpace,
         DataTransferKit::HostSpace>( success, out );
}
#endif

//---------------------------------------------------------------------------//
// end tstMapInterface.cpp
//---------------------------------------------------------------------------//
