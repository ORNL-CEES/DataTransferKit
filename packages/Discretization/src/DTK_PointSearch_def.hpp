/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_POINT_SEARCH_DEF_HPP
#define DTK_POINT_SEARCH_DEF_HPP

#include <DTK_DBC.hpp>
#include <DTK_Topology.hpp>

#include <Teuchos_SerializationTraits.hpp>

namespace Teuchos
{
template <typename Ordinal>
class SerializationTraits<Ordinal, DataTransferKit::Point>
    : public DirectSerializationTraits<Ordinal, DataTransferKit::Point>
{
};
}

namespace DataTransferKit
{
namespace internal
{
template <typename DeviceType>
KOKKOS_FUNCTION void computeBlockCellsBoundingBox(
    unsigned int const dim, int const i, unsigned int const n_nodes,
    unsigned int const node_offset, unsigned int const topo_id,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<unsigned int *, DeviceType> offset,
    Kokkos::View<double **, DeviceType> coordinates,
    Kokkos::View<double ***, DeviceType> block_cells,
    Kokkos::View<Box *, DeviceType> bounding_boxes,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell )
{
    Box bounding_box;
    // If dim == 2, we need to set bounding_box.minCorner()[2] and
    // bounding_box.maxCorner[2].
    if ( dim == 2 )
    {
        bounding_box.minCorner()[2] = 0;
        bounding_box.maxCorner()[2] = 1;
    }
    unsigned int const k = offset( i );
    for ( unsigned int node = 0; node < n_nodes; ++node )
    {
        unsigned int const n = node_offset + node;
        for ( unsigned int d = 0; d < dim; ++d )
        {
            // Copy the coordinated in block_cells
            block_cells( k, node, d ) = coordinates( cells( n ), d );
            // Build the bounding box.
            if ( block_cells( k, node, d ) < bounding_box.minCorner()[d] )
                bounding_box.minCorner()[d] = block_cells( k, node, d );
            if ( block_cells( k, node, d ) > bounding_box.maxCorner()[d] )
                bounding_box.maxCorner()[d] = block_cells( k, node, d );
        }
    }
    bounding_boxes( i ) = bounding_box;
    bounding_box_to_cell( i, topo_id ) = k;
}
}

template <typename DeviceType>
PointSearch<DeviceType>::PointSearch(
    Teuchos::RCP<const Teuchos::Comm<int>> comm,
    Kokkos::View<DTK_CellTopology *, DeviceType> const &cell_topologies_view,
    Kokkos::View<unsigned int *, DeviceType> const &cells,
    Kokkos::View<double **, DeviceType> const &nodes_coordinates,
    Kokkos::View<double **, DeviceType> const &points_coordinates,
    std::string const &strategy )
    : _comm( comm )
    , _target_to_source_distributor( _comm )
    , _bounding_boxes( "bounding_boxes", cell_topologies_view.extent( 0 ) )
    , _bounding_box_to_cell( "bounding_box_to_cell",
                             cell_topologies_view.extent( 0 ), DTK_N_TOPO )
{
    // Initialize _bounding_box_to_cell to an invalid state
    Kokkos::deep_copy( _bounding_box_to_cell, static_cast<unsigned int>( -1 ) );

    // Compute the number of cells of each of the supported topologies.
    std::array<unsigned int, DTK_N_TOPO> n_cells_per_topo =
        computeTopologies( cell_topologies_view );

    // Convert the cells and nodes_coordinates View to block_cells
    convertMesh( n_cells_per_topo, cell_topologies_view, cells,
                 nodes_coordinates );

    // Perform the distributed search
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> per_topo_ranks;
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<Point *, DeviceType> imported_points( "imported_points", 0 );
    Kokkos::View<int *, DeviceType> imported_query_ids( "imported_query_ids",
                                                        0 );
    Kokkos::View<int *, DeviceType> imported_cell_indices( "imported_indices",
                                                           0 );
    Kokkos::View<int *, DeviceType> ranks( "ranks", 0 );
    // The distributed search only supports dim == 3, so we need to convert the
    // 2D points to 3D points.
    if ( _dim == 2 )
    {
        unsigned int const n_points = points_coordinates.extent( 0 );
        Kokkos::View<double **, DeviceType> points_coord_3d( "points_coord_3d",
                                                             n_points, 3 );
        Kokkos::parallel_for(
            REGION_NAME( "convert_2D_pts_to_3D_pts" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
            KOKKOS_LAMBDA( int const i ) {
                points_coord_3d( i, 0 ) = points_coordinates( i, 0 );
                points_coord_3d( i, 1 ) = points_coordinates( i, 1 );
                points_coord_3d( i, 2 ) = 0.;
            } );
        Kokkos::fence();

        performDistributedSearch( points_coord_3d, imported_points,
                                  imported_query_ids, imported_cell_indices,
                                  ranks );
    }
    else
        performDistributedSearch( points_coordinates, imported_points,
                                  imported_query_ids, imported_cell_indices,
                                  ranks );

    // We need to separate the data for the different topologies because of
    // Intrepid2. Because a point can be found in multiple cells, we need to
    // compute the number of cells of each topology.
    unsigned int const n_imports = imported_points.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> topo( "topo", n_imports );
    Kokkos::View<unsigned int *, DeviceType> topo_size( "topo_size",
                                                        DTK_N_TOPO );
    // We want to use _bounding_box_to_cell in the lambda function but we
    // can't because it is a private variable so create a local variable
    // instead.
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell =
        _bounding_box_to_cell;
    Kokkos::parallel_for(
        "separate_topo", Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
        KOKKOS_LAMBDA( int const i ) {
            for ( unsigned int j = 0; j < DTK_N_TOPO; ++j )
            {
                if ( bounding_box_to_cell( imported_cell_indices( i ), j ) !=
                     static_cast<unsigned int>( -1 ) )
                {
                    topo( i ) = j;
                    Kokkos::atomic_increment( &topo_size( j ) );
                    break;
                }
            }
        } );
    Kokkos::fence();

    // Now that we know the size, allocate all the View.
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_cell_indices;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_query_ids;
    std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_reference_points;
    std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_point_in_cell;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_ranks;
    auto topo_size_host = Kokkos::create_mirror_view( topo_size );
    Kokkos::deep_copy( topo_size_host, topo_size );
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        filtered_per_topo_cell_indices[i] = Kokkos::View<int *, DeviceType>(
            "filtered_per_topo_cell_indices_" + std::to_string( i ),
            topo_size_host( i ) );
        filtered_per_topo_query_ids[i] = Kokkos::View<int *, DeviceType>(
            "filtered_per_topo_query_ids_" + std::to_string( i ),
            topo_size_host( i ) );
        filtered_per_topo_reference_points[i] =
            Kokkos::View<double **, DeviceType>(
                "filtered_per_topo_reference_points_" + std::to_string( i ),
                topo_size_host( i ), _dim );
        filtered_per_topo_point_in_cell[i] = Kokkos::View<bool *, DeviceType>(
            "filtered_per_topo_point_in_cell_" + std::to_string( i ),
            topo_size_host( i ) );
        filtered_per_topo_ranks[i] = Kokkos::View<int *, DeviceType>(
            "filtered_per_topo_ranks_" + std::to_string( i ),
            topo_size_host( i ) );
    }

    // Check if the points are in the cells
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        if ( _block_cells[topo_id].extent( 0 ) != 0 )
        {
            Kokkos::View<double **, DeviceType> filtered_points(
                "filtered_points", topo_size_host( topo_id ), _dim );
            Kokkos::View<unsigned int *, DeviceType> point_indices_map(
                "point_indices_map", topo_size_host( topo_id ) );
            performPointInCell(
                _block_cells[topo_id], imported_cell_indices, imported_points,
                imported_query_ids, ranks, topo, topo_id, filtered_points,
                filtered_per_topo_cell_indices[topo_id],
                filtered_per_topo_query_ids[topo_id], point_indices_map,
                filtered_per_topo_reference_points[topo_id],
                filtered_per_topo_point_in_cell[topo_id],
                filtered_per_topo_ranks[topo_id] );
        }

    // Filter the points
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> filtered_ranks;
    filter( strategy, filtered_per_topo_point_in_cell,
            filtered_per_topo_reference_points, filtered_per_topo_cell_indices,
            filtered_per_topo_query_ids, filtered_per_topo_ranks,
            filtered_ranks );

    // Build the _source_to_target_distributor
    build_distributor( filtered_ranks );
}

template <typename DeviceType>
std::tuple<Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
           Kokkos::View<Point *, DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>>
PointSearch<DeviceType>::getSearchResults()
{
    // Flatten the results
    unsigned int size = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        size += _reference_points[topo_id].extent( 0 );

    Kokkos::View<int *, DeviceType> ranks( "ranks", size );
    Kokkos::View<int *, DeviceType> cell_indices( "cell_indices", size );
    Kokkos::View<Point *, DeviceType> ref_pts( "ref_pts", size );
    Kokkos::View<unsigned int *, DeviceType> query_ids( "query_ids", size );

    // FIXME to do on the device
    Kokkos::deep_copy( ranks, _comm->getRank() );
    auto cell_indices_host = Kokkos::create_mirror_view( cell_indices );
    Kokkos::deep_copy( cell_indices_host, cell_indices );
    auto ref_pts_host = Kokkos::create_mirror_view( ref_pts );
    Kokkos::deep_copy( ref_pts_host, ref_pts );
    auto query_ids_host = Kokkos::create_mirror_view( query_ids );
    Kokkos::deep_copy( query_ids_host, query_ids );

    unsigned int pos = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int const size = _query_ids[topo_id].extent( 0 );
        for ( unsigned int i = 0; i < size; ++i )
        {
            cell_indices_host( pos ) = _cell_indices[topo_id]( i );
            query_ids_host( pos ) = _query_ids[topo_id]( i );
            for ( unsigned int d = 0; d < _dim; ++d )
                ref_pts_host( pos )[d] = _reference_points[topo_id]( i, d );
            ++pos;
        }
    }
    Kokkos::deep_copy( cell_indices, cell_indices_host );
    Kokkos::deep_copy( query_ids, query_ids_host );

    // Communicate the results
    unsigned int n_imports =
        _target_to_source_distributor.getTotalReceiveLength();
    Kokkos::View<int *, DeviceType> imported_ranks( "imported_ranks",
                                                    n_imports );
    Kokkos::View<int *, DeviceType> imported_cell_indices(
        "imported_cell_indices", n_imports );
    Kokkos::View<Point *, DeviceType> imported_ref_pts( "imported_ref_pts",
                                                        n_imports );
    Kokkos::View<unsigned int *, DeviceType> imported_query_ids(
        "imported_query_ids", n_imports );

    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, ranks, imported_ranks );
    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, cell_indices, imported_cell_indices );
    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, ref_pts, imported_ref_pts );
    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, query_ids, imported_query_ids );

    Kokkos::View<unsigned int *, DeviceType> offset_view( "offset_view", 0 );
    if ( n_imports != 0 )
    {
        // Sort the results by query ids
        // Because of the MPI communications and the sorting by topologies, all
        // the query have been reordered. So we put them back in the initial
        // order using the query ids.
        using ExecutionSpace = typename DeviceType::execution_space;
        typedef Kokkos::BinOp1D<Kokkos::View<unsigned int *, DeviceType>>
            CompType;
        Kokkos::Experimental::MinMaxScalar<unsigned int> result;
        Kokkos::Experimental::MinMax<unsigned int> reducer( result );
        parallel_reduce(
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
            Kokkos::Impl::min_max_functor<
                Kokkos::View<unsigned int *, DeviceType>>( imported_query_ids ),
            reducer );
        if ( result.min_val != result.max_val )
        {
            Kokkos::BinSort<Kokkos::View<unsigned int *, DeviceType>, CompType>
                bin_sort(
                    imported_query_ids,
                    CompType( n_imports / 2, result.min_val, result.max_val ),
                    true );
            bin_sort.create_permute_vector();
            bin_sort.sort( imported_query_ids );
            bin_sort.sort( imported_cell_indices );
            bin_sort.sort( imported_ranks );
            bin_sort.sort( imported_ref_pts );
        }

        // Compute the offset
        // FIXME to do on the device
        auto imported_query_ids_host =
            Kokkos::create_mirror_view( imported_query_ids );
        std::vector<unsigned int> offset;
        offset.push_back( 0 );
        for ( unsigned int i = 1; i < imported_query_ids.extent( 0 ); ++i )
        {
            if ( imported_query_ids_host( i ) !=
                 imported_query_ids_host( i - 1 ) )
                offset.push_back( i );
        }
        offset.push_back( imported_query_ids.extent( 0 ) );
        Kokkos::realloc( offset_view, offset.size() );
        auto offset_view_host = Kokkos::create_mirror_view( offset_view );
        for ( unsigned int i = 0; i < offset.size(); ++i )
            offset_view_host( i ) = offset[i];
        Kokkos::deep_copy( offset_view, offset_view_host );
    }
    else
    {
        Kokkos::realloc( offset_view, 1 );
        Kokkos::deep_copy( offset_view, 0 );
    }

    return std::make_tuple( imported_ranks, imported_cell_indices,
                            imported_ref_pts, offset_view );
}

template <typename DeviceType>
std::array<unsigned int, DTK_N_TOPO> PointSearch<DeviceType>::computeTopologies(
    Kokkos::View<DTK_CellTopology *, DeviceType> const &cell_topologies_view )
{
    // This needs to be done on the host because we will use n_cells_per_topo
    // to allocate Kokkos::View and this needs to be done on the host
    std::array<unsigned int, DTK_N_TOPO> n_cells_per_topo;
    auto cell_topologies_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_host, cell_topologies_view );
    unsigned int const n_cells = cell_topologies_view.extent( 0 );
    DTK_CellTopology dtk_cell_topo;
    n_cells_per_topo.fill( 0 );
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        dtk_cell_topo = cell_topologies_host( i );
        n_cells_per_topo[dtk_cell_topo] += 1;
    }

    // We do not support meshes that contains both 2D and 3D cells. All the
    // cells are either 2D or 3D
    Topologies topologies;
    _dim = topologies[dtk_cell_topo].dim;
#if HAVE_DTK_DBC
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        if ( n_cells_per_topo[i] != 0 )
            DTK_REQUIRE( topologies[dtk_cell_topo].dim == _dim );
    }
#endif

    // PointInCell does not support all Intrepid2 topologies
    DTK_REQUIRE( n_cells_per_topo[DTK_TET_11] == 0 );
    DTK_REQUIRE( n_cells_per_topo[DTK_HEX_20] == 0 );
    DTK_REQUIRE( n_cells_per_topo[DTK_WEDGE_15] == 0 );

    return n_cells_per_topo;
}

template <typename DeviceType>
void PointSearch<DeviceType>::convertMesh(
    std::array<unsigned int, DTK_N_TOPO> const &n_cells_per_topo,
    Kokkos::View<DTK_CellTopology *, DeviceType> const &cell_topologies_view,
    Kokkos::View<unsigned int *, DeviceType> const &cells,
    Kokkos::View<double **, DeviceType> const &coordinates )
{
    // First, we get the number of nodes for each topology to get initialize
    // _block_cells at the right size.
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo(
        "n_nodes_per_topo" );
    auto n_nodes_per_topo_host = Kokkos::create_mirror_view( n_nodes_per_topo );
    Topologies topologies;
    for ( int i = 0; i < DTK_N_TOPO; ++i )
        n_nodes_per_topo_host( i ) = topologies[i].n_nodes;
    Kokkos::deep_copy( n_nodes_per_topo, n_nodes_per_topo_host );

    // Create an array of Kokkos::View where each View is a block of cells with
    // the same topology
    for ( int i = 0; i < DTK_N_TOPO; ++i )
    {
        _block_cells[i] = Kokkos::View<double ***, DeviceType>(
            "block_cells_" + std::to_string( i ), n_cells_per_topo[i],
            n_nodes_per_topo_host( i ), _dim );
    }

    // Compute the offset associated to each cell in the coordinates View.
    unsigned int const n_cells = cell_topologies_view.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> nodes_per_cell( "nodes_per_cell",
                                                             n_cells );
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view =
        cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> node_offset( "node_offset",
                                                          n_cells );
    computeNodeOffset( n_cells, cell_topo_view, n_nodes_per_topo,
                       nodes_per_cell, node_offset );

    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<unsigned int *, DeviceType> offset( "offset", n_cells );
    unsigned int const dim = _dim;
    Kokkos::View<Box *, DeviceType> bounding_boxes = _bounding_boxes;
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell =
        _bounding_box_to_cell;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        // Compute the relative position of each cell of a given topology
        computeOffset( cell_topo_view, topo_id, offset );

        // Fill _block_cell and create the bounding boxes from _block_cells
        Kokkos::View<double ***, DeviceType> block_cells =
            _block_cells[topo_id];
        Kokkos::parallel_for(
            REGION_NAME( "fill_block_cell_" + std::to_string( topo_id ) ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
            KOKKOS_LAMBDA( int const i ) {
                if ( cell_topo_view( i ) == topo_id )
                {
                    internal::computeBlockCellsBoundingBox(
                        dim, i, n_nodes_per_topo( topo_id ), node_offset( i ),
                        topo_id, cells, offset, coordinates, block_cells,
                        bounding_boxes, bounding_box_to_cell );
                }
            } );
        Kokkos::fence();
    }
}

template <typename DeviceType>
void PointSearch<DeviceType>::performDistributedSearch(
    Kokkos::View<double **, DeviceType> points_coord,
    Kokkos::View<Point *, DeviceType> &imported_points,
    Kokkos::View<int *, DeviceType> &imported_query_ids,
    Kokkos::View<int *, DeviceType> &imported_cell_indices,
    Kokkos::View<int *, DeviceType> &ranks )
{
    DistributedSearchTree<DeviceType> distributed_tree( _comm,
                                                        _bounding_boxes );

    unsigned int const n_points = points_coord.extent( 0 );

    // Build the queries
    // FIXME do not use Within predicate because it requires a radius, i.e., it
    // is mesh dependent
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<Details::Within *, DeviceType> queries( "queries", n_points );
    Kokkos::parallel_for(
        "register_queries", Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
        KOKKOS_LAMBDA( int i ) {
            queries( i ) =
                Details::within( {{points_coord( i, 0 ), points_coord( i, 1 ),
                                   points_coord( i, 2 )}},
                                 1e-14 );
        } );
    Kokkos::fence();

    // Perform the distributed search
    Kokkos::View<int *, DeviceType> indices( "indices" );
    Kokkos::View<int *, DeviceType> offset( "offset" );
    distributed_tree.query( queries, indices, offset, ranks );

    // Create the source to target distributor
    auto ranks_host = Kokkos::create_mirror_view( ranks );
    Kokkos::deep_copy( ranks_host, ranks );
    Tpetra::Distributor source_to_target_distributor( _comm );
    unsigned int const n_imports = source_to_target_distributor.createFromSends(
        Teuchos::ArrayView<int const>( ranks_host.data(),
                                       ranks_host.extent( 0 ) ) );

    // Communicate cell indices
    Kokkos::realloc( imported_cell_indices, n_imports );
    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, indices, imported_cell_indices );
    // Duplicate the points_coord for the communication. Duplicating the points
    // allows us to use the same distributor.
    unsigned int const indices_size = indices.extent( 0 );
    Kokkos::View<Point *, DeviceType> exported_points( "exported_points",
                                                       indices_size );
    Kokkos::View<int *, DeviceType> exported_query_ids( "exported_query_ids",
                                                        indices_size );
    // This line should not be necessary but there is problem with the
    // lambda capture on CUDA otherwise.
    unsigned int dim = _dim;
    Kokkos::parallel_for(
        "duplicate_points",
        Kokkos::RangePolicy<ExecutionSpace>( 0, offset.extent( 0 ) - 1 ),
        KOKKOS_LAMBDA( int const i ) {
            for ( int j = offset( i ); j < offset( i + 1 ); ++j )
            {
                exported_query_ids( j ) = i;
                for ( unsigned int k = 0; k < dim; ++k )
                    exported_points( j )[k] = points_coord( i, k );
            }
        } );
    Kokkos::fence();

    // Communicate the points
    Kokkos::realloc( imported_points, n_imports );
    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, exported_points, imported_points );

    // Communicate the query_ids. We communicate the query_ids to keep track of
    // which points is associated to which query.
    Kokkos::realloc( imported_query_ids, n_imports );
    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, exported_query_ids, imported_query_ids );

    // Communicate the ranks of the sending processors. This will be used to
    // build the _target_to_source_distributor.
    Kokkos::realloc( ranks, n_imports );
    Kokkos::View<int *, DeviceType> exported_ranks( "exported_ranks",
                                                    indices_size );
    Kokkos::deep_copy( exported_ranks, _comm->getRank() );
    DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, exported_ranks, ranks );
}

template <typename DeviceType>
void PointSearch<DeviceType>::filterTopology(
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    Kokkos::View<int *, DeviceType> cell_indices,
    Kokkos::View<Point *, DeviceType> points,
    Kokkos::View<int *, DeviceType> query_ids,
    Kokkos::View<int *, DeviceType> ranks,
    Kokkos::View<unsigned int *, DeviceType> point_indices_map,
    Kokkos::View<int *, DeviceType> filtered_cell_indices,
    Kokkos::View<double **, DeviceType> filtered_points,
    Kokkos::View<int *, DeviceType> filtered_query_ids,
    Kokkos::View<int *, DeviceType> filtered_ranks )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const n_imports = topo.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> offset( "offset", n_imports );
    computeOffset( topo, topo_id, offset );

    // Create Kokkos::View with the points and the cell indices associated
    // with cells of topo_id topology. Also transform 3D points back to 2D
    // points.
    // We want to use _bounding_box_to_cell in the lambda function but we want
    // can't because it is a private variable so create a local variable
    // instead.
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell =
        _bounding_box_to_cell;
    unsigned int dim = _dim;
    Kokkos::parallel_for(
        "filter_data", Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
        KOKKOS_LAMBDA( int const i ) {
            if ( topo( i ) == topo_id )
            {
                unsigned int const k = offset( i );
                filtered_cell_indices( k ) =
                    bounding_box_to_cell( cell_indices( i ), topo_id );
                for ( unsigned int j = 0; j < dim; ++j )
                    filtered_points( k, j ) = points( i )[j];
                filtered_query_ids( k ) = query_ids( i );
                point_indices_map( k ) = i;
                filtered_ranks( k ) = ranks( i );
            }
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void PointSearch<DeviceType>::computeNodeOffset(
    unsigned int const n_cells,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view,
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
    Kokkos::View<unsigned int *, DeviceType> nodes_per_cell,
    Kokkos::View<unsigned int *, DeviceType> node_offset )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for( REGION_NAME( "fill_nodes_per_cell" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
                          KOKKOS_LAMBDA( int const i ) {
                              nodes_per_cell( i ) =
                                  n_nodes_per_topo( cell_topo_view( i ) );
                          } );
    Kokkos::fence();
    Kokkos::parallel_scan(
        REGION_NAME( "compute_node_offset" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
        KOKKOS_LAMBDA( int i, int &update, bool final_pass ) {
            if ( final_pass )
                node_offset( i ) = update;
            update += nodes_per_cell( i );
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void PointSearch<DeviceType>::performPointInCell(
    Kokkos::View<double ***, DeviceType> cells,
    Kokkos::View<int *, DeviceType> imported_cell_indices,
    Kokkos::View<Point *, DeviceType> imported_points,
    Kokkos::View<int *, DeviceType> imported_query_ids,
    Kokkos::View<int *, DeviceType> imported_ranks,
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    Kokkos::View<double **, DeviceType> filtered_points,
    Kokkos::View<int *, DeviceType> filtered_cell_indices,
    Kokkos::View<int *, DeviceType> filtered_query_ids,
    Kokkos::View<unsigned int *, DeviceType> point_indices_map,
    Kokkos::View<double **, DeviceType> reference_points,
    Kokkos::View<bool *, DeviceType> point_in_cell,
    Kokkos::View<int *, DeviceType> filtered_ranks )
{
    filterTopology( topo, topo_id, imported_cell_indices, imported_points,
                    imported_query_ids, imported_ranks, point_indices_map,
                    filtered_cell_indices, filtered_points, filtered_query_ids,
                    filtered_ranks );

    // Perform the PointInCell search
    Topologies topologies;
    PointInCell<DeviceType>::search(
        filtered_points, cells, filtered_cell_indices, topologies[topo_id].topo,
        reference_points, point_in_cell );
}

template <typename DeviceType>
void PointSearch<DeviceType>::filter(
    std::string const &strategy,
    std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_point_in_cell,
    std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_reference_points,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_cell_indices,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_query_ids,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const &ranks,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> &filtered_ranks )
{
    filter_all( filtered_per_topo_point_in_cell,
                filtered_per_topo_reference_points,
                filtered_per_topo_cell_indices, filtered_per_topo_query_ids,
                ranks, filtered_ranks );
    //    if ( strategy == "unique" )
    //        // Need to do an extra step to get ride of the duplicated point
    //        filter_unique();
}

template <typename DeviceType>
void PointSearch<DeviceType>::filter_all(
    std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_point_in_cell,
    std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_reference_points,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_cell_indices,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_query_ids,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const &ranks,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> &filtered_ranks )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int dim = _dim;

    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int n_ref_points =
            filtered_per_topo_point_in_cell[topo_id].extent( 0 );
        if ( n_ref_points != 0 )
        {
            int n_filtered_ref_points = 0;
            Kokkos::View<bool *, DeviceType> pt_in_cell =
                filtered_per_topo_point_in_cell[topo_id];
            Kokkos::parallel_reduce(
                REGION_NAME( "compute_n_ref_pts" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
                KOKKOS_LAMBDA( int i, int &partial_sum ) {
                    if ( pt_in_cell[i] == true )
                        partial_sum += 1;
                },
                n_filtered_ref_points );

            // We are only interested in points that belong to the cells. So we
            // need to filter out all the points that were false positive of
            // the distributed search.
            Kokkos::realloc( _reference_points[topo_id], n_filtered_ref_points,
                             _dim );
            Kokkos::realloc( _query_ids[topo_id], n_filtered_ref_points );
            Kokkos::realloc( _cell_indices[topo_id], n_filtered_ref_points );
            Kokkos::realloc( filtered_ranks[topo_id], n_filtered_ref_points );

            // We cannot use private member in a lambda function with CUDA
            Kokkos::View<Coordinate **, DeviceType> ref_points =
                _reference_points[topo_id];
            Kokkos::View<int *, DeviceType> query_ids = _query_ids[topo_id];
            Kokkos::View<int *, DeviceType> cell_indices =
                _cell_indices[topo_id];

            Kokkos::View<unsigned int *, DeviceType> offset( "offset",
                                                             n_ref_points );
            computeOffset( pt_in_cell, true, offset );
            Kokkos::parallel_for(
                REGION_NAME( "filter" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
                KOKKOS_LAMBDA( int const i ) {
                    if ( pt_in_cell[i] )
                    {
                        unsigned int k = offset( i );
                        for ( unsigned int d = 0; d < dim; ++d )
                            ref_points( k, d ) =
                                filtered_per_topo_reference_points[topo_id](
                                    i, d );
                        query_ids( k ) =
                            filtered_per_topo_query_ids[topo_id]( i );
                        cell_indices( k ) =
                            filtered_per_topo_cell_indices[topo_id]( i );
                        filtered_ranks[topo_id]( k ) = ranks[topo_id]( i );
                    }
                } );
            Kokkos::fence();
        }
    }
}

template <typename DeviceType>
void PointSearch<DeviceType>::build_distributor(
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_ranks )
{
    // Flatten the filtered ranks to be used by the distributor
    std::vector<int> flatten_ranks;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        auto rank_host = Kokkos::create_mirror_view( filtered_ranks[topo_id] );
        Kokkos::deep_copy( rank_host, filtered_ranks[topo_id] );
        unsigned int const rank_host_size = rank_host.size();
        for ( unsigned int i = 0; i < rank_host_size; ++i )
            flatten_ranks.push_back( rank_host( i ) );
    }

    _target_to_source_distributor.createFromSends(
        Teuchos::ArrayView<int const>( flatten_ranks ) );
}
}

// Explicit instantiation macro
#define DTK_POINTSEARCH_INSTANT( NODE )                                        \
    template class PointSearch<typename NODE::device_type>;

#endif
