/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include "DTK_ConfigDefs.hpp"
#include "DTK_Types.h"

#include "PointCloudProblemGenerator/ExodusProblemGenerator.hpp"
#include "PointCloudProblemGenerator/PointCloudProblemGenerator.hpp"
#include <ArborX.hpp>
#include <DTK_NearestNeighborOperator.hpp>
#include <DTK_ParallelTraits.hpp>

#include <Kokkos_View.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_Export.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <cmath>
#include <limits>
#include <string>

//---------------------------------------------------------------------------//
// Brute-force calculate the nearest neighbors
template <class Device, class... CoordViewProperties>
void computeNeighborsBruteForce(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
    const Kokkos::View<DataTransferKit::Coordinate **, CoordViewProperties...>
        &src_coords,
    const Kokkos::View<DataTransferKit::Coordinate **, CoordViewProperties...>
        &tgt_coords,
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, Device>
        &src_nearest_coords )
{
    // Copy coordinates to the host.
    size_t num_src = src_coords.extent( 0 );
    size_t num_tgt = tgt_coords.extent( 0 );
    size_t space_dim = src_coords.extent( 1 );
    auto host_sources = Kokkos::create_mirror_view( src_coords );
    Kokkos::deep_copy( host_sources, src_coords );
    auto host_targets = Kokkos::create_mirror_view( tgt_coords );
    Kokkos::deep_copy( host_targets, tgt_coords );

    // Get the number of source points on each rank.
    int comm_size = comm->getSize();
    int comm_rank = comm->getRank();
    Teuchos::Array<size_t> rank_num_src( comm_size, 0 );
    rank_num_src[comm_rank] = num_src;
    Teuchos::Array<size_t> local_sizes( comm_size, 0 );
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, comm_size,
                        rank_num_src.getRawPtr(), local_sizes.getRawPtr() );

    // Compute the global number of source points.
    size_t num_global = 0;
    for ( const auto &l : local_sizes )
        num_global += l;

    // Calculate rank offsets.
    Teuchos::Array<size_t> offsets( comm_size, 0 );
    for ( int i = 1; i < comm_size; ++i )
    {
        offsets[i] = offsets[i - 1] + local_sizes[i - 1];
    }

    // Allocate a send and receive view for creating a global list of source
    // coordinates.
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, Device>
        rank_src_coords( "rank_src_coords", num_global, space_dim );
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, Device>
        all_src_coords( "src_nearest_coords", num_global, space_dim );

    // Copy the local rank source coordinates into the appropriate spot in the
    // send buffer for the global source coordinate view.
    for ( size_t d = 0; d < space_dim; ++d )
    {
        auto src_subview = Kokkos::subview( host_sources, Kokkos::ALL, d );
        auto send_subview = Kokkos::subview( rank_src_coords, Kokkos::ALL, d );
        std::copy( src_subview.data(), src_subview.data() + num_src,
                   send_subview.data() + offsets[comm_rank] );
    }

    // Reduce source coordinate data across all ranks to build a
    // globally-replicated source coordinate view.
    for ( size_t d = 0; d < space_dim; ++d )
    {
        auto send_subview = Kokkos::subview( rank_src_coords, Kokkos::ALL, d );
        auto recv_subview = Kokkos::subview( all_src_coords, Kokkos::ALL, d );
        Teuchos::reduceAll<int, DataTransferKit::Coordinate>(
            *comm, Teuchos::REDUCE_SUM, Teuchos::as<int>( num_global ),
            send_subview.data(), recv_subview.data() );
    }

    // Allocate the nearest neighbors array.
    src_nearest_coords =
        Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft,
                     Device>( "src_nearest_coords", num_tgt, space_dim );

    // Use a brute-force method to find the nearest neighbors. Compute the
    // distance for each target one at a time compared to all source points.
    double distance;
    for ( size_t t = 0; t < num_tgt; ++t )
    {
        // Reset the minimum distance for this target.
        double min_distance = std::numeric_limits<double>::max();

        // Check each source to see if it is the closest.
        for ( size_t s = 0; s < num_global; ++s )
        {
            // Compute the distance between the source and target.
            distance = 0.0;
            for ( size_t d = 0; d < space_dim; ++d )
                distance += ( host_targets( t, d ) - all_src_coords( s, d ) ) *
                            ( host_targets( t, d ) - all_src_coords( s, d ) );
            distance = std::sqrt( distance );

            // If the distance is a new minimum for this target set the
            // source coordinates as the new nearest and update the new
            // minimum distance.
            if ( distance < min_distance )
            {
                for ( size_t d = 0; d < space_dim; ++d )
                    src_nearest_coords( t, d ) = all_src_coords( s, d );
                min_distance = distance;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Test a problem that is uniquely owned. Returns a view of ints in the
// target decomposition indicating success/failure for each point.
template <class... ViewProperties>
void testUniquelyOwnedProblem(
    const Kokkos::View<DataTransferKit::Coordinate **, ViewProperties...>
        &src_coords,
    const Kokkos::View<double **, ViewProperties...> &src_field,
    const Kokkos::View<DataTransferKit::Coordinate **, ViewProperties...>
        &tgt_coords,
    const Kokkos::View<double **, ViewProperties...> &tgt_field, bool &success,
    Teuchos::FancyOStream &out )
{
    // NOTE following line is there to get rid of unused parameter warning.
    // However, I am not sure why this is a parameter at all since the exact
    // answer is computed via brute force search of the nearest neighbor and
    // applying the field function space onto the coordinates of that nearest
    // neighbor.
    std::ignore = tgt_field;

    // Types.
    using CoordView =
        Kokkos::View<DataTransferKit::Coordinate **, ViewProperties...>;
    using Device = typename CoordView::device_type;
    using ExecutionSpace = typename Device::execution_space;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Sizes
    size_t space_dim = src_coords.extent( 1 );
    size_t num_src = src_coords.extent( 0 );
    size_t num_tgt = tgt_coords.extent( 0 );

    // Data function.
    auto data_func = KOKKOS_LAMBDA( const double x, const double y,
                                    const double z, const int comp )
    {
        return x * ( comp + 1 ) + y * ( comp + 2 ) + z * ( comp + 3 );
    };

    // Generate a source field.
    auto fill_src = KOKKOS_LAMBDA( const size_t i )
    {
        src_field( i, 0 ) = data_func( src_coords( i, 0 ), src_coords( i, 1 ),
                                       src_coords( i, 2 ), 0 );
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, num_src ),
                          fill_src );
    Kokkos::fence();

    // For now we copy the coordinates to the device default layout for
    // compatibility with the current operator definition.
    Kokkos::View<DataTransferKit::Coordinate **, Device> src_coords_device(
        "src_coords_device", num_src, space_dim );
    Kokkos::deep_copy( src_coords_device, src_coords );

    Kokkos::View<DataTransferKit::Coordinate **, Device> tgt_coords_device(
        "tgt_coords_device", num_tgt, space_dim );
    Kokkos::deep_copy( tgt_coords_device, tgt_coords );

    // Create a nearest neighbor operator.
    DataTransferKit::NearestNeighborOperator<Device> nearest_op(
        *( Teuchos::rcp_dynamic_cast<Teuchos::MpiComm<int> const>( comm )
               ->getRawMpiComm() ),
        src_coords_device, tgt_coords_device );

    // For now we copy the source field and create a target field with views
    // that use the default layout for the device type. Also note that we are
    // only using a single field dimension to be compatible with the current
    // operator definition.
    Kokkos::View<double *, Device> src_field_device( "src_field_device",
                                                     num_src );
    auto src_field_copy = KOKKOS_LAMBDA( const size_t i )
    {
        src_field_device( i ) = src_field( i, 0 );
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, num_src ),
                          src_field_copy );
    Kokkos::fence();

    Kokkos::View<double *, Device> tgt_field_device( "tgt_field_device",
                                                     num_tgt );

    // Apply the operator.
    nearest_op.apply( src_field_device, tgt_field_device );

    // Create the expected nearest coordinates.
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, Device>
        nearest_src_coords;
    computeNeighborsBruteForce( comm, src_coords, tgt_coords,
                                nearest_src_coords );

    // Check that the field evaluated at the points found via brute-force
    // search is the same as those computed by the operator.
    auto host_tgt_field = Kokkos::create_mirror_view( tgt_field_device );
    Kokkos::deep_copy( host_tgt_field, tgt_field_device );
    double data_val;
    for ( size_t t = 0; t < num_tgt; ++t )
    {
        data_val =
            data_func( nearest_src_coords( t, 0 ), nearest_src_coords( t, 1 ),
                       nearest_src_coords( t, 2 ), 0 );
        TEST_FLOATING_EQUALITY( host_tgt_field( t ), data_val, 1.0e-12 );
    }
}

//---------------------------------------------------------------------------//
// Partition the grids one-to-one
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ExodusProblemGenerator, one_to_one, Node )
{
    // Type aliases.
    using DeviceType = typename Node::device_type;
    using Scalar = double;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Create a problem generator. Two tet meshes of a sphere centered at
    // (0,0,0) with a radius of 10 are used. The source mesh is represented by
    // a fine tet mesh and the target mesh is represented by a coarse tet mesh.
    std::string src_exodus_file = "fine_sphere.exo";
    std::string tgt_exodus_file = "coarse_sphere.exo";
    DataTransferKit::ExodusProblemGenerator<Scalar, DeviceType, DeviceType>
        generator(
            *( Teuchos::rcp_dynamic_cast<Teuchos::MpiComm<int> const>( comm )
                   ->getRawMpiComm() ),
            src_exodus_file, tgt_exodus_file );

    // Generate a uniquely owned problem.
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        src_coords;
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        tgt_coords;
    Kokkos::View<Scalar **, Kokkos::LayoutLeft, DeviceType> src_field;
    Kokkos::View<Scalar **, Kokkos::LayoutLeft, DeviceType> tgt_field;
    generator.createUniquelyOwnedProblem( src_coords, src_field, tgt_coords,
                                          tgt_field );

    // Test the problem.
    testUniquelyOwnedProblem( src_coords, src_field, tgt_coords, tgt_field,
                              success, out );
}

//---------------------------------------------------------------------------//
// Partition the grids with standard ghosting from connectivity.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ExodusProblemGenerator, ghosted, Node )
{
    // Type aliases.
    using DeviceType = typename Node::device_type;
    using Scalar = double;
    using TpetraMap = Tpetra::Map<int, DataTransferKit::GlobalOrdinal, Node>;
    using TpetraCoordVector =
        Tpetra::MultiVector<DataTransferKit::Coordinate, int,
                            DataTransferKit::GlobalOrdinal, Node>;
    using TpetraScalarVector =
        Tpetra::MultiVector<Scalar, int, DataTransferKit::GlobalOrdinal, Node>;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Create a problem generator. Two tet meshes of a sphere centered at
    // (0,0,0) with a radius of 10 are used. The source mesh is represented by
    // a fine tet mesh and the target mesh is represented by a coarse tet mesh.
    std::string src_exodus_file = "fine_sphere.exo";
    std::string tgt_exodus_file = "coarse_sphere.exo";
    DataTransferKit::ExodusProblemGenerator<Scalar, DeviceType, DeviceType>
        generator(
            *( Teuchos::rcp_dynamic_cast<Teuchos::MpiComm<int> const>( comm )
                   ->getRawMpiComm() ),
            src_exodus_file, tgt_exodus_file );

    // Generate a ghosted owned problem.
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        src_coords;
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        tgt_coords;
    Kokkos::View<DataTransferKit::GlobalOrdinal *, Kokkos::LayoutLeft,
                 DeviceType>
        src_gids;
    Kokkos::View<DataTransferKit::GlobalOrdinal *, Kokkos::LayoutLeft,
                 DeviceType>
        tgt_gids;
    Kokkos::View<Scalar **, Kokkos::LayoutLeft, DeviceType> src_field;
    Kokkos::View<Scalar **, Kokkos::LayoutLeft, DeviceType> tgt_field;
    generator.createGhostedProblem( src_coords, src_gids, src_field, tgt_coords,
                                    tgt_gids, tgt_field );

    // Create Tpetra maps from the problem global ids.
    auto invalid_ordinal =
        Teuchos::OrdinalTraits<DataTransferKit::GlobalOrdinal>::invalid();
    Teuchos::RCP<const TpetraMap> src_ghost_map =
        Teuchos::rcp( new TpetraMap( invalid_ordinal, src_gids, 0, comm ) );
    Teuchos::RCP<const TpetraMap> tgt_ghost_map =
        Teuchos::rcp( new TpetraMap( invalid_ordinal, tgt_gids, 0, comm ) );

    // Create Tpetra vectors from the ghosted problem coordinates.
    TpetraCoordVector src_coords_ghost( src_ghost_map, src_coords );
    TpetraCoordVector tgt_coords_ghost( tgt_ghost_map, tgt_coords );
    TpetraScalarVector src_field_ghost( src_ghost_map, src_field );
    TpetraScalarVector tgt_field_ghost( tgt_ghost_map, tgt_field );

    // Create unique Tpetra vectors for coordinates and fields.
    auto src_unique_map = Tpetra::createOneToOne( src_ghost_map );
    auto tgt_unique_map = Tpetra::createOneToOne( tgt_ghost_map );
    int space_dim = src_coords.extent( 1 );
    int field_ncomp = src_field.extent( 1 );
    TpetraCoordVector src_coords_unique( src_unique_map, space_dim );
    TpetraCoordVector tgt_coords_unique( tgt_unique_map, space_dim );
    TpetraScalarVector src_field_unique( src_unique_map, field_ncomp );
    TpetraScalarVector tgt_field_unique( tgt_unique_map, field_ncomp );

    // Export the ghosted problem coordinates to a uniquely owned
    // decomposition.
    auto src_ghost_to_unique_exporter =
        Tpetra::createExport( src_ghost_map, src_unique_map );
    auto tgt_ghost_to_unique_exporter =
        Tpetra::createExport( tgt_ghost_map, tgt_unique_map );
    src_coords_unique.doExport( src_coords_ghost, *src_ghost_to_unique_exporter,
                                Tpetra::REPLACE );
    tgt_coords_unique.doExport( tgt_coords_ghost, *tgt_ghost_to_unique_exporter,
                                Tpetra::REPLACE );
    src_field_unique.doExport( src_field_ghost, *src_ghost_to_unique_exporter,
                               Tpetra::REPLACE );
    tgt_field_unique.doExport( tgt_field_ghost, *tgt_ghost_to_unique_exporter,
                               Tpetra::REPLACE );

    // Test the uniquely owned problem.
    auto src_coords_unique_view =
        src_coords_unique.template getLocalView<DeviceType>();
    auto tgt_coords_unique_view =
        tgt_coords_unique.template getLocalView<DeviceType>();
    auto src_field_unique_view =
        src_field_unique.template getLocalView<DeviceType>();
    auto tgt_field_unique_view =
        tgt_field_unique.template getLocalView<DeviceType>();
    testUniquelyOwnedProblem( src_coords_unique_view, src_field_unique_view,
                              tgt_coords_unique_view, tgt_field_unique_view,
                              success, out );
}

//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ExodusProblemGenerator, one_to_one,  \
                                          NODE )                               \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ExodusProblemGenerator, ghosted,     \
                                          NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
