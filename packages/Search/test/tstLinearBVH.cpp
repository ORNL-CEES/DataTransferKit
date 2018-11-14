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

#include <DTK_LinearBVH.hpp>

#include <Teuchos_UnitTestHarness.hpp>

#include "DTK_BoostRTreeHelpers.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <tuple>

#include "Search_UnitTestHelpers.hpp"

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, empty_tree, DeviceType )
{
    // tree is empty, it has no leaves.
    for ( auto const &empty_bvh : {
              DataTransferKit::BVH<DeviceType>{}, // default constructed
              makeBvh<DeviceType>( {} ), // constructed with empty view of boxes
          } )
    {
        TEST_ASSERT( empty_bvh.empty() );
        TEST_EQUALITY( empty_bvh.size(), 0 );
        // BVH::bounds() returns an invalid box when the tree is empty.
        TEST_ASSERT(
            DataTransferKit::Details::equals( empty_bvh.bounds(), {} ) );

        // Passing a view with no query does seem a bit silly but we still need
        // to support it. And since the tag dispatching yields different tree
        // traversals for nearest and spatial predicates, we do have to check
        // the results for various type of queries.
        checkResults( empty_bvh, makeOverlapQueries<DeviceType>( {} ), {}, {0},
                      success, out );

        // NOTE: Admittedly testing for both overlap and within queries might be
        // a bit overkill but I'd rather test for all the queries we plan on
        // using.
        checkResults( empty_bvh, makeWithinQueries<DeviceType>( {} ), {}, {0},
                      success, out );

        checkResults( empty_bvh, makeNearestQueries<DeviceType>( {} ), {}, {0},
                      success, out );

        // Passing an empty distance vector.
        checkResults( empty_bvh, makeNearestQueries<DeviceType>( {} ), {}, {0},
                      {}, success, out );

        // Now passing a couple queries of various type and checking the
        // results.
        checkResults(
            empty_bvh,
            makeOverlapQueries<DeviceType>( {
                {}, // Did not bother giving a valid box here but that's fine.
                {},
            } ),
            {}, {0, 0, 0}, success, out );

        checkResults( empty_bvh,
                      makeWithinQueries<DeviceType>( {
                          {{{0., 0., 0.}}, 1.},
                          {{{1., 1., 1.}}, 2.},
                      } ),
                      {}, {0, 0, 0}, success, out );

        checkResults( empty_bvh,
                      makeNearestQueries<DeviceType>( {
                          {{{0., 0., 0.}}, 1},
                          {{{1., 1., 1.}}, 2},
                      } ),
                      {}, {0, 0, 0}, success, out );

        checkResults( empty_bvh,
                      makeNearestQueries<DeviceType>( {
                          {{{0., 0., 0.}}, 1},
                          {{{1., 1., 1.}}, 2},
                      } ),
                      {}, {0, 0, 0}, {}, success, out );
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, single_leaf_tree, DeviceType )
{
    // tree has a single leaf (unit box)
    auto const bvh = makeBvh<DeviceType>( {
        {{{0., 0., 0.}}, {{1., 1., 1.}}},
    } );

    TEST_ASSERT( !bvh.empty() );
    TEST_EQUALITY( bvh.size(), 1 );
    TEST_ASSERT( DataTransferKit::Details::equals(
        bvh.bounds(), {{{0., 0., 0.}}, {{1., 1., 1.}}} ) );

    checkResults( bvh, makeOverlapQueries<DeviceType>( {} ), {}, {0}, success,
                  out );

    checkResults( bvh, makeWithinQueries<DeviceType>( {} ), {}, {0}, success,
                  out );

    checkResults( bvh, makeNearestQueries<DeviceType>( {} ), {}, {0}, success,
                  out );

    checkResults( bvh, makeNearestQueries<DeviceType>( {} ), {}, {0}, {},
                  success, out );

    checkResults( bvh,
                  makeOverlapQueries<DeviceType>( {
                      {{{5., 5., 5.}}, {{5., 5., 5.}}},
                      {{{.5, .5, .5}}, {{.5, .5, .5}}},
                  } ),
                  {0}, {0, 0, 1}, success, out );

    checkResults( bvh,
                  makeWithinQueries<DeviceType>( {
                      {{{0., 0., 0.}}, 1.},
                      {{{1., 1., 1.}}, 3.},
                      {{{5., 5., 5.}}, 2.},
                  } ),
                  {0, 0}, {0, 1, 2, 2}, success, out );

    checkResults( bvh,
                  makeNearestQueries<DeviceType>( {
                      {{{0., 0., 0.}}, 1},
                      {{{1., 1., 1.}}, 2},
                      {{{2., 2., 2.}}, 3},
                  } ),
                  {0, 0, 0}, {0, 1, 2, 3}, success, out );

    checkResults( bvh,
                  makeNearestQueries<DeviceType>( {
                      {{{1., 0., 0.}}, 1},
                      {{{0., 2., 0.}}, 2},
                      {{{0., 0., 3.}}, 3},
                  } ),
                  {0, 0, 0}, {0, 1, 2, 3}, {0., 1., 2.}, success, out );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, couple_leaves_tree, DeviceType )
{
    auto const bvh = makeBvh<DeviceType>( {
        {{{0., 0., 0.}}, {{0., 0., 0.}}},
        {{{1., 1., 1.}}, {{1., 1., 1.}}},
    } );

    TEST_ASSERT( !bvh.empty() );
    TEST_EQUALITY( bvh.size(), 2 );
    TEST_ASSERT( DataTransferKit::Details::equals(
        bvh.bounds(), {{{0., 0., 0.}}, {{1., 1., 1.}}} ) );

    // single query overlap with nothing
    checkResults( bvh,
                  makeOverlapQueries<DeviceType>( {
                      {},
                  } ),
                  {}, {0, 0}, success, out );

    // single query overlap with both
    checkResults( bvh,
                  makeOverlapQueries<DeviceType>( {
                      {{{0., 0., 0.}}, {{1., 1., 1.}}},
                  } ),
                  {1, 0}, {0, 2}, success, out );

    // single query overlap with only one
    checkResults( bvh,
                  makeOverlapQueries<DeviceType>( {
                      {{{0.5, 0.5, 0.5}}, {{1.5, 1.5, 1.5}}},
                  } ),
                  {1}, {0, 1}, success, out );

    // a couple queries both overlap with nothing
    checkResults( bvh,
                  makeOverlapQueries<DeviceType>( {
                      {},
                      {},
                  } ),
                  {}, {0, 0, 0}, success, out );

    // a couple queries first overlap with nothing second with only one
    checkResults( bvh,
                  makeOverlapQueries<DeviceType>( {
                      {},
                      {{{0., 0., 0.}}, {{0., 0., 0.}}},
                  } ),
                  {0}, {0, 0, 1}, success, out );

    // no query
    checkResults( bvh, makeOverlapQueries<DeviceType>( {} ), {}, {0}, success,
                  out );

    checkResults( bvh,
                  makeNearestQueries<DeviceType>( {
                      {{{0., 0., 0.}}, 2},
                      {{{1., 0., 0.}}, 4},
                  } ),
                  {0, 1, 0, 1}, {0, 2, 4}, {0., sqrt( 3. ), 1., sqrt( 2. )},
                  success, out );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, duplicated_leaves, DeviceType )
{
    // The tree contains multiple (more than two) leaves that will be assigned
    // the same Morton code.  This was able to trigger a bug that we discovered
    // when building trees over ~10M indexable values.  The hierarchy generated
    // at construction had leaves with no parent which yielded a segfault later
    // when computing bounding boxes and walking the hierarchy toward the root.
    auto const bvh = makeBvh<DeviceType>( {
        {{{0., 0., 0.}}, {{0., 0., 0.}}},
        {{{1., 1., 1.}}, {{1., 1., 1.}}},
        {{{1., 1., 1.}}, {{1., 1., 1.}}},
        {{{1., 1., 1.}}, {{1., 1., 1.}}},
    } );

    // I commented out the test below because it will fail unless we modify
    // checkResults() such that it sorts the results for each query before
    // comparing arrays because there is no way to predict the ordering of the
    // indices for the duplicate leaves.  This is fine, the point of this test
    // was to fix the segfault we would get when calling the constructor of BVH.
    // checkResults( bvh,
    //              makeWithinQueries<DeviceType>( {
    //                  {{{0., 0., 0.}}, 1.},
    //                  {{{1., 1., 1.}}, 1.},
    //                  {{{.5, .5, .5}}, 1.},
    //              } ),
    //              {0, 1, 2, 3, 0, 1, 2, 3}, {0, 1, 4, 8}, success, out );

    // Suppress a warning about success and out not being used
    TEST_EQUALITY( true, true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, buffer_optimization, DeviceType )
{
    auto const bvh = makeBvh<DeviceType>( {
        {{{0., 0., 0.}}, {{0., 0., 0.}}},
        {{{1., 0., 0.}}, {{1., 0., 0.}}},
        {{{2., 0., 0.}}, {{2., 0., 0.}}},
        {{{3., 0., 0.}}, {{3., 0., 0.}}},
    } );

    auto const queries = makeOverlapQueries<DeviceType>( {
        {},
        {{{0., 0., 0.}}, {{3., 3., 3.}}},
        {},
    } );

    using ViewType = Kokkos::View<int *, DeviceType>;
    ViewType indices( "indices" );
    ViewType offset( "offset" );

    std::vector<int> indices_ref = {3, 2, 1, 0};
    std::vector<int> offset_ref = {0, 0, 4, 4};
    auto checkResultsAreFine = [&indices, &offset, &indices_ref, &offset_ref,
                                &success, &out]() -> void {
        TEST_COMPARE_ARRAYS( indices, indices_ref );
        TEST_COMPARE_ARRAYS( offset, offset_ref );
    };

    TEST_NOTHROW( bvh.query( queries, indices, offset ) );
    checkResultsAreFine();

    // compute number of results per query
    auto counts = DataTransferKit::cloneWithoutInitializingNorCopying( offset );
    DataTransferKit::adjacentDifference( offset, counts );
    // extract optimal buffer size
    auto const max_results_per_query = DataTransferKit::max( counts );
    TEST_EQUALITY( max_results_per_query, 4 );

    // optimal size
    TEST_NOTHROW(
        bvh.query( queries, indices, offset, -max_results_per_query ) );
    checkResultsAreFine();

    // buffer size insufficient
    TEST_COMPARE( max_results_per_query, >, 1 );
    TEST_NOTHROW( bvh.query( queries, indices, offset, +1 ) );
    checkResultsAreFine();
    TEST_THROW( bvh.query( queries, indices, offset, -1 ),
                DataTransferKit::DataTransferKitException );

    // adequate buffer size
    TEST_COMPARE( max_results_per_query, <, 5 );
    TEST_NOTHROW( bvh.query( queries, indices, offset, +5 ) );
    checkResultsAreFine();
    TEST_NOTHROW( bvh.query( queries, indices, offset, -5 ) );
    checkResultsAreFine();

    // passing null size skips the buffer optimization and never throws
    TEST_NOTHROW( bvh.query( queries, indices, offset, 0 ) );
    checkResultsAreFine();
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, not_exceeding_stack_capacity,
                                   DeviceType )
{
    std::vector<DataTransferKit::Box> boxes;
    int const n = 4096; // exceed stack capacity which is 64
    boxes.reserve( n );
    for ( int i = 0; i < n; ++i )
    {
        double const a = i;
        double const b = i + 1;
        boxes.push_back( {{{a, a, a}}, {{b, b, b}}} );
    }
    auto const bvh = makeBvh<DeviceType>( boxes );

    Kokkos::View<int *, DeviceType> indices( "indices" );
    Kokkos::View<int *, DeviceType> offset( "offset" );
    // query number of nearest neighbors that exceed capacity of the stack is
    // not a problem
    TEST_NOTHROW( bvh.query( makeNearestQueries<DeviceType>( {
                                 {{{0., 0., 0.}}, n},
                             } ),
                             indices, offset ) );
    TEST_EQUALITY( DataTransferKit::lastElement( offset ), n );

    // spatial query that find all indexable in the tree is also fine
    TEST_NOTHROW( bvh.query( makeOverlapQueries<DeviceType>( {
                                 {},
                                 {{{0., 0., 0.}}, {{n, n, n}}},
                             } ),
                             indices, offset ) );
    TEST_EQUALITY( DataTransferKit::lastElement( offset ), n );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, miscellaneous, DeviceType )
{
    auto const bvh = makeBvh<DeviceType>( {
        {{{1., 3., 5.}}, {{2., 4., 6.}}},
    } );
    auto const empty_bvh = makeBvh<DeviceType>( {} );

    // Batched queries BVH::query( Kokkos::View<Query *, ...>, ... ) returns
    // early if the tree is empty.  Below we ensure that a direct call to the
    // single query TreeTraversal::query() actually handles empty trees
    // properly.
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<int *, DeviceType> zeros( "zeros", 3 );
    Kokkos::deep_copy( zeros, 255 );
    Kokkos::View<Kokkos::pair<int, double> *, DeviceType> empty_buffer(
        "empty_buffer" );
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecutionSpace>( 0, 1 ), KOKKOS_LAMBDA( int ) {
            DataTransferKit::Point p = {{0., 0., 0.}};
            double r = 1.0;
            // spatial query on empty tree
            zeros( 0 ) =
                DataTransferKit::Details::TreeTraversal<DeviceType>::query(
                    empty_bvh, DataTransferKit::within( p, r ), []( int ) {} );
            // nearest query on empty tree
            zeros( 1 ) =
                DataTransferKit::Details::TreeTraversal<DeviceType>::query(
                    empty_bvh, DataTransferKit::nearest( p ),
                    []( int, double ) {}, empty_buffer );
            // nearest query for k < 1
            zeros( 2 ) =
                DataTransferKit::Details::TreeTraversal<DeviceType>::query(
                    bvh, DataTransferKit::nearest( p, 0 ), []( int, double ) {},
                    empty_buffer );
        } );
    Kokkos::fence();
    auto zeros_host = Kokkos::create_mirror_view( zeros );
    Kokkos::deep_copy( zeros_host, zeros );
    std::vector<int> zeros_ref = {0, 0, 0};
    TEST_COMPARE_ARRAYS( zeros_host, zeros_ref );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, structured_grid, DeviceType )
{
    double Lx = 100.0;
    double Ly = 100.0;
    double Lz = 100.0;
    int nx = 11;
    int ny = 11;
    int nz = 11;
    int n = nx * ny * nz;

    using ExecutionSpace = typename DeviceType::execution_space;

    Kokkos::View<DataTransferKit::Box *, DeviceType> bounding_boxes(
        "bounding_boxes", n );
    Kokkos::parallel_for(
        "fill_bounding_boxes", Kokkos::RangePolicy<ExecutionSpace>( 0, nx ),
        KOKKOS_LAMBDA( int i ) {
            [[gnu::unused]] double x, y, z;
            for ( int j = 0; j < ny; ++j )
                for ( int k = 0; k < nz; ++k )
                {
                    x = i * Lx / ( nx - 1 );
                    y = j * Ly / ( ny - 1 );
                    z = k * Lz / ( nz - 1 );
                    bounding_boxes[i + j * nx + k * ( nx * ny )] = {
                        {{x, y, z}}, {{x, y, z}}};
                }
        } );
    Kokkos::fence();

    DataTransferKit::BVH<DeviceType> bvh( bounding_boxes );

    // (i) use same objects for the queries than the objects we constructed the
    // BVH
    // i-2  i-1  i  i+1
    //
    //  o    o   o   o   j+1
    //          ---
    //  o    o | x | o   j
    //          ---
    //  o    o   o   o   j-1
    //
    //  o    o   o   o   j-2
    //
    Kokkos::View<DataTransferKit::Overlap *, DeviceType> queries( "queries",
                                                                  n );
    Kokkos::parallel_for(
        "fill_queries", Kokkos::RangePolicy<ExecutionSpace>( 0, n ),
        KOKKOS_LAMBDA( int i ) {
            queries( i ) = DataTransferKit::Overlap( bounding_boxes( i ) );
        } );
    Kokkos::fence();

    Kokkos::View<int *, DeviceType> indices( "indices", n );
    Kokkos::View<int *, DeviceType> offset( "offset", n );

    bvh.query( queries, indices, offset );

    auto indices_host = Kokkos::create_mirror_view( indices );
    Kokkos::deep_copy( indices_host, indices );
    auto offset_host = Kokkos::create_mirror_view( offset );
    Kokkos::deep_copy( offset_host, offset );

    // we expect the collision list to be diag(0, 1, ..., nx*ny*nz-1)
    for ( int i = 0; i < n; ++i )
    {
        TEST_EQUALITY( indices_host( i ), i );
        TEST_EQUALITY( offset_host( i ), i );
    }

    // (ii) use bounding boxes that overlap with first neighbors
    //
    // i-2  i-1  i  i+1
    //
    //  o    x---x---x   j+1
    //       |       |
    //  o    x   x   x   j
    //       |       |
    //  o    x---x---x   j-1
    //
    //  o    o   o   o   j-2
    //

    auto bounding_boxes_host = Kokkos::create_mirror_view( bounding_boxes );
    std::function<int( int, int, int )> ind = [nx, ny]( int i, int j, int k ) {
        return i + j * nx + k * ( nx * ny );
    };
    std::vector<std::set<int>> ref( n );
    for ( int i = 0; i < nx; ++i )
        for ( int j = 0; j < ny; ++j )
            for ( int k = 0; k < nz; ++k )
            {
                int const index = ind( i, j, k );
                // bounding box around nodes of the structured grid will
                // overlap with neighboring nodes
                bounding_boxes_host[index] = {
                    {{( i - 1 ) * Lx / ( nx - 1 ), ( j - 1 ) * Ly / ( ny - 1 ),
                      ( k - 1 ) * Lz / ( nz - 1 )}},
                    {{( i + 1 ) * Lx / ( nx - 1 ), ( j + 1 ) * Ly / ( ny - 1 ),
                      ( k + 1 ) * Lz / ( nz - 1 )}}};
                // fill in reference solution to check againt the collision
                // list computed during the tree traversal
                if ( ( i > 0 ) && ( j > 0 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i - 1, j - 1, k - 1 ) );
                if ( ( i > 0 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i - 1, j, k - 1 ) );
                if ( ( i > 0 ) && ( j < ny - 1 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i - 1, j + 1, k - 1 ) );
                if ( ( i > 0 ) && ( j > 0 ) )
                    ref[index].emplace( ind( i - 1, j - 1, k ) );
                if ( i > 0 )
                    ref[index].emplace( ind( i - 1, j, k ) );
                if ( ( i > 0 ) && ( j < ny - 1 ) )
                    ref[index].emplace( ind( i - 1, j + 1, k ) );
                if ( ( i > 0 ) && ( j > 0 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i - 1, j - 1, k + 1 ) );
                if ( ( i > 0 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i - 1, j, k + 1 ) );
                if ( ( i > 0 ) && ( j < ny - 1 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i - 1, j + 1, k + 1 ) );

                if ( ( j > 0 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i, j - 1, k - 1 ) );
                if ( k > 0 )
                    ref[index].emplace( ind( i, j, k - 1 ) );
                if ( ( j < ny - 1 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i, j + 1, k - 1 ) );
                if ( j > 0 )
                    ref[index].emplace( ind( i, j - 1, k ) );
                if ( true )
                    ref[index].emplace( ind( i, j, k ) );
                if ( j < ny - 1 )
                    ref[index].emplace( ind( i, j + 1, k ) );
                if ( ( j > 0 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i, j - 1, k + 1 ) );
                if ( k < nz - 1 )
                    ref[index].emplace( ind( i, j, k + 1 ) );
                if ( ( j < ny - 1 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i, j + 1, k + 1 ) );

                if ( ( i < nx - 1 ) && ( j > 0 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i + 1, j - 1, k - 1 ) );
                if ( ( i < nx - 1 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i + 1, j, k - 1 ) );
                if ( ( i < nx - 1 ) && ( j < ny - 1 ) && ( k > 0 ) )
                    ref[index].emplace( ind( i + 1, j + 1, k - 1 ) );
                if ( ( i < nx - 1 ) && ( j > 0 ) )
                    ref[index].emplace( ind( i + 1, j - 1, k ) );
                if ( i < nx - 1 )
                    ref[index].emplace( ind( i + 1, j, k ) );
                if ( ( i < nx - 1 ) && ( j < ny - 1 ) )
                    ref[index].emplace( ind( i + 1, j + 1, k ) );
                if ( ( i < nx - 1 ) && ( j > 0 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i + 1, j - 1, k + 1 ) );
                if ( ( i < nx - 1 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i + 1, j, k + 1 ) );
                if ( ( i < nx - 1 ) && ( j < ny - 1 ) && ( k < nz - 1 ) )
                    ref[index].emplace( ind( i + 1, j + 1, k + 1 ) );
            }

    Kokkos::deep_copy( bounding_boxes, bounding_boxes_host );
    Kokkos::parallel_for(
        "fill_first_neighbors_queries",
        Kokkos::RangePolicy<ExecutionSpace>( 0, n ), KOKKOS_LAMBDA( int i ) {
            queries[i] = DataTransferKit::Overlap( bounding_boxes[i] );
        } );
    Kokkos::fence();
    bvh.query( queries, indices, offset );
    Kokkos::deep_copy( indices_host, indices );
    Kokkos::deep_copy( offset_host, offset );

    for ( int i = 0; i < nx; ++i )
        for ( int j = 0; j < ny; ++j )
            for ( int k = 0; k < nz; ++k )
            {
                int index = ind( i, j, k );
                for ( int l = offset( index ); l < offset( index + 1 ); ++l )
                {
                    TEST_ASSERT( ref[index].count( indices( l ) ) != 0 );
                }
            }

    // (iii) use random points
    //
    // i-1      i      i+1
    //
    //  o       o       o   j+1
    //         -------
    //        |       |
    //        |   +   |
    //  o     | x     | o   j
    //         -------
    //
    //  o       o       o   j-1
    //
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution_x( 0.0, Lz );
    std::uniform_real_distribution<double> distribution_y( 0.0, Ly );
    std::uniform_real_distribution<double> distribution_z( 0.0, Lz );

    for ( int l = 0; l < n; ++l )
    {
        double x = distribution_x( generator );
        double y = distribution_y( generator );
        double z = distribution_z( generator );
        bounding_boxes_host( l ) = {
            {{x - 0.5 * Lx / ( nx - 1 ), y - 0.5 * Ly / ( ny - 1 ),
              z - 0.5 * Lz / ( nz - 1 )}},
            {{x + 0.5 * Lx / ( nx - 1 ), y + 0.5 * Ly / ( ny - 1 ),
              z + 0.5 * Lz / ( nz - 1 )}}};

        int i = std::round( x / Lx * ( nx - 1 ) );
        int j = std::round( y / Ly * ( ny - 1 ) );
        int k = std::round( z / Lz * ( nz - 1 ) );
        // Save the indices for the check
        ref[l] = {ind( i, j, k )};
    }

    Kokkos::deep_copy( bounding_boxes, bounding_boxes_host );
    Kokkos::parallel_for(
        "fill_first_neighbors_queries",
        Kokkos::RangePolicy<ExecutionSpace>( 0, n ), KOKKOS_LAMBDA( int i ) {
            queries[i] = DataTransferKit::Overlap( bounding_boxes[i] );
        } );
    Kokkos::fence();
    bvh.query( queries, indices, offset );
    Kokkos::deep_copy( indices_host, indices );
    Kokkos::deep_copy( offset_host, offset );

    for ( int i = 0; i < n; ++i )
    {
        TEST_EQUALITY( offset( i ), i );
        TEST_ASSERT( ref[i].count( indices[i] ) != 0 );
    }
}

std::vector<std::array<double, 3>>
make_stuctured_cloud( double Lx, double Ly, double Lz, int nx, int ny, int nz )
{
    std::vector<std::array<double, 3>> cloud( nx * ny * nz );
    std::function<int( int, int, int )> ind = [nx, ny]( int i, int j, int k ) {
        return i + j * nx + k * ( nx * ny );
    };
    double x, y, z;
    for ( int i = 0; i < nx; ++i )
        for ( int j = 0; j < ny; ++j )
            for ( int k = 0; k < nz; ++k )
            {
                x = i * Lx / ( nx - 1 );
                y = j * Ly / ( ny - 1 );
                z = k * Lz / ( nz - 1 );
                cloud[ind( i, j, k )] = {{x, y, z}};
            }
    return cloud;
}

std::vector<std::array<double, 3>> make_random_cloud( double Lx, double Ly,
                                                      double Lz, int n )
{
    std::vector<std::array<double, 3>> cloud( n );
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution_x( 0.0, Lx );
    std::uniform_real_distribution<double> distribution_y( 0.0, Ly );
    std::uniform_real_distribution<double> distribution_z( 0.0, Lz );
    for ( int i = 0; i < n; ++i )
    {
        double x = distribution_x( generator );
        double y = distribution_y( generator );
        double z = distribution_z( generator );
        cloud[i] = {{x, y, z}};
    }
    return cloud;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LinearBVH, rtree, DeviceType )
{
    // contruct a cloud of points (nodes of a structured grid)
    double Lx = 10.0;
    double Ly = 10.0;
    double Lz = 10.0;
    int nx = 11;
    int ny = 11;
    int nz = 11;
    auto cloud = make_stuctured_cloud( Lx, Ly, Lz, nx, ny, nz );
    int n = cloud.size();

    Kokkos::View<DataTransferKit::Box *, DeviceType> bounding_boxes(
        "bounding_boxes", n );
    auto bounding_boxes_host = Kokkos::create_mirror_view( bounding_boxes );
    // build bounding volume hierarchy
    for ( int i = 0; i < n; ++i )
    {
        auto const &point = cloud[i];
        double x = std::get<0>( point );
        double y = std::get<1>( point );
        double z = std::get<2>( point );
        bounding_boxes_host[i] = {{{x, y, z}}, {{x, y, z}}};
    }

    Kokkos::deep_copy( bounding_boxes, bounding_boxes_host );

    DataTransferKit::BVH<DeviceType> bvh( bounding_boxes );

    auto rtree = BoostRTreeHelpers::makeRTree( bounding_boxes );

    // random points for radius search and kNN queries
    // compare our solution against Boost R-tree
    int const n_points = 100;
    auto queries = make_random_cloud( Lx, Ly, Lz, n_points );
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<double * [3], ExecutionSpace> point_coords( "point_coords",
                                                             n_points );
    auto point_coords_host = Kokkos::create_mirror_view( point_coords );
    Kokkos::View<double *, ExecutionSpace> radii( "radii", n_points );
    auto radii_host = Kokkos::create_mirror_view( radii );
    Kokkos::View<int * [2], ExecutionSpace> within_n_pts( "within_n_pts",
                                                          n_points );
    Kokkos::View<int * [2], ExecutionSpace> nearest_n_pts( "nearest_n_pts",
                                                           n_points );
    Kokkos::View<int *, ExecutionSpace> k( "distribution_k", n_points );
    auto k_host = Kokkos::create_mirror_view( k );
    // use random radius for the search and random number k of for the kNN
    // search
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution_radius(
        0.0, std::sqrt( Lx * Lx + Ly * Ly + Lz * Lz ) );
    std::uniform_int_distribution<int> distribution_k(
        1, std::floor( sqrt( nx * nx + ny * ny + nz * nz ) ) );
    for ( unsigned int i = 0; i < n_points; ++i )
    {
        auto const &point = queries[i];
        double x = std::get<0>( point );
        double y = std::get<1>( point );
        double z = std::get<2>( point );
        radii_host[i] = distribution_radius( generator );
        k_host[i] = distribution_k( generator );
        point_coords_host( i, 0 ) = x;
        point_coords_host( i, 1 ) = y;
        point_coords_host( i, 2 ) = z;
    }

    Kokkos::deep_copy( point_coords, point_coords_host );
    Kokkos::deep_copy( radii, radii_host );
    Kokkos::deep_copy( k, k_host );

    Kokkos::View<DataTransferKit::Nearest<DataTransferKit::Point> *, DeviceType>
        nearest_queries( "nearest_queries", n_points );
    Kokkos::parallel_for(
        "register_nearest_queries",
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
        KOKKOS_LAMBDA( int i ) {
            nearest_queries( i ) =
                DataTransferKit::nearest<DataTransferKit::Point>(
                    {{point_coords( i, 0 ), point_coords( i, 1 ),
                      point_coords( i, 2 )}},
                    k( i ) );
        } );
    Kokkos::fence();

    Kokkos::View<DataTransferKit::Within *, DeviceType> within_queries(
        "within_queries", n_points );
    Kokkos::parallel_for( "register_within_queries",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
                          KOKKOS_LAMBDA( int i ) {
                              within_queries( i ) = DataTransferKit::within(
                                  {{point_coords( i, 0 ), point_coords( i, 1 ),
                                    point_coords( i, 2 )}},
                                  radii( i ) );
                          } );
    Kokkos::fence();

    auto rtree_results =
        BoostRTreeHelpers::performQueries( rtree, nearest_queries );

    Kokkos::View<int *, DeviceType> offset_nearest( "offset_nearest" );
    Kokkos::View<int *, DeviceType> indices_nearest( "indices_nearest" );
    bvh.query( nearest_queries, indices_nearest, offset_nearest );
    auto bvh_results = std::make_tuple( offset_nearest, indices_nearest );

    validateResults( rtree_results, bvh_results, success, out );

    Kokkos::View<int *, DeviceType> offset_within( "offset_within" );
    Kokkos::View<int *, DeviceType> indices_within( "indices_within" );
    bvh.query( within_queries, indices_within, offset_within );
    bvh_results = std::make_tuple( offset_within, indices_within );

    rtree_results = BoostRTreeHelpers::performQueries( rtree, within_queries );

    validateResults( rtree_results, bvh_results, success, out );
}

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, empty_tree,               \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, single_leaf_tree,         \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, couple_leaves_tree,       \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, duplicated_leaves,        \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, buffer_optimization,      \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        LinearBVH, not_exceeding_stack_capacity, DeviceType##NODE )            \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, miscellaneous,            \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, structured_grid,          \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LinearBVH, rtree, DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
