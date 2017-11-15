/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_INTERPOLATION_DEF_HPP
#define DTK_INTERPOLATION_DEF_HPP

#include <DTK_InterpolationUtils.hpp>
#include <DTK_PointInCell.hpp>

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
template <typename DeviceType>
Interpolation<DeviceType>::Interpolation(
    Teuchos::RCP<const Teuchos::Comm<int>> comm,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<double **, DeviceType> nodes_coordinates,
    Kokkos::View<double **, DeviceType> points_coordinates, DTK_FEType fe_type )
    : _apply_allowed( false )
    , _comm( comm )
    , _distributor( comm )
    , _cell_topologies_view( cell_topologies_view )
    , _bounding_boxes( "bounding_boxes", _cell_topologies_view.extent( 0 ) )
    , _bounding_box_to_cell( "bounding_box_to_cell",
                             _cell_topologies_view.extent( 0 ), DTK_N_TOPO )
{
    DTK_INSIST( fe_type == DTK_HGRAD );

    _cell_topologies[0] = shards::getCellTopologyData<shards::Triangle<3>>();
    _cell_topologies[1] = shards::getCellTopologyData<shards::Triangle<6>>();
    _cell_topologies[2] =
        shards::getCellTopologyData<shards::Quadrilateral<4>>();
    _cell_topologies[3] =
        shards::getCellTopologyData<shards::Quadrilateral<9>>();
    _cell_topologies[4] = shards::getCellTopologyData<shards::Tetrahedron<4>>();
    _cell_topologies[5] =
        shards::getCellTopologyData<shards::Tetrahedron<10>>();
    _cell_topologies[6] =
        shards::getCellTopologyData<shards::Tetrahedron<11>>();
    _cell_topologies[7] = shards::getCellTopologyData<shards::Hexahedron<8>>();
    _cell_topologies[8] = shards::getCellTopologyData<shards::Hexahedron<20>>();
    _cell_topologies[9] = shards::getCellTopologyData<shards::Hexahedron<27>>();
    _cell_topologies[10] = shards::getCellTopologyData<shards::Pyramid<5>>();
    _cell_topologies[11] = shards::getCellTopologyData<shards::Pyramid<13>>();
    _cell_topologies[12] = shards::getCellTopologyData<shards::Wedge<6>>();
    _cell_topologies[13] = shards::getCellTopologyData<shards::Wedge<15>>();
    _cell_topologies[14] = shards::getCellTopologyData<shards::Wedge<18>>();

    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        _n_cells_per_topo[topo_id] = 0;

    findReferencePoints( cells, nodes_coordinates, points_coordinates );
}

template <typename DeviceType>
Interpolation<DeviceType>::Interpolation(
    Teuchos::RCP<const Teuchos::Comm<int>> comm,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<double **, DeviceType> nodes_coordinates,
    Kokkos::View<double **, DeviceType> points_coordinates,
    Kokkos::View<LocalOrdinal *, DeviceType> object_dof_ids,
    DTK_FEType fe_type )
    : Interpolation<DeviceType>( comm, cell_topologies_view, cells,
                                 nodes_coordinates, points_coordinates,
                                 fe_type )
{
    _apply_allowed = true;

    // This function does exactly the same thing as convertMesh but for
    // object_dof_ids

    // Number of nodes for each topology.
    // shards objects don't exist on the GPU, so we need to work on the CPU and
    // then copy the data back
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo(
        "n_nodes_per_topo" );
    auto n_nodes_per_topo_host = Kokkos::create_mirror_view( n_nodes_per_topo );
    for ( int i = 0; i < DTK_N_TOPO; ++i )
        n_nodes_per_topo_host( i ) = _cell_topologies[i].getNodeCount();
    Kokkos::deep_copy( n_nodes_per_topo, n_nodes_per_topo_host );

    // Create an array of Kokkos::View where each View is a block of cells with
    // the same topology
    for ( int i = 0; i < DTK_N_TOPO; ++i )
    {
        _block_cells[i] = Kokkos::View<double ***, DeviceType>(
            "block_cells_" + std::to_string( i ), _n_cells_per_topo[i],
            n_nodes_per_topo_host( i ), _dim );
    }

    // Fill _block_cells from the cells and coordinates
    unsigned int const n_cells = _cell_topologies_view.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> nodes_per_cell( "nodes_per_cell",
                                                             n_cells );
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view =
        _cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> offset( "offset", n_cells );
    Kokkos::View<unsigned int *, DeviceType> node_offset( "node_offset",
                                                          n_cells );
    computeNodeOffset( n_cells, cell_topo_view, n_nodes_per_topo,
                       nodes_per_cell, node_offset );

    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int const n_local_cells = _n_cells_per_topo[topo_id];
        if ( n_local_cells != 0 )
        {
            unsigned int const n_nodes_per_topo =
                _cell_topologies[topo_id].getNodeCount();
            // Assume that the basis functions match the cell topologies
            _cell_dofs_ids[topo_id] = Kokkos::View<LocalOrdinal **, DeviceType>(
                "cell_dofs_ids_" + std::to_string( topo_id ), n_local_cells,
                n_nodes_per_topo );

            Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids =
                _cell_dofs_ids[topo_id];

            computeOffset( cell_topologies_view, topo_id, offset );

            fillCellDofIds( topo_id, n_cells, n_nodes_per_topo,
                            cell_topologies_view, offset, node_offset,
                            object_dof_ids, cell_dofs_ids );
        }
    }
}

template <typename DeviceType>
Interpolation<DeviceType>::Interpolation(
    Teuchos::RCP<const Teuchos::Comm<int>> comm,
    Kokkos::View<DTK_CellTopology *, DeviceType>,
    Kokkos::View<unsigned int *, DeviceType>,
    Kokkos::View<double **, DeviceType>, Kokkos::View<double **, DeviceType>,
    Kokkos::View<unsigned int *, DeviceType>,
    Kokkos::View<DTK_Quadrature *, DeviceType>, DTK_FEType )
    : _distributor( comm )
{
    DTK_INSIST( false );
}

template <typename DeviceType>
Interpolation<DeviceType>::Interpolation(
    Teuchos::RCP<const Teuchos::Comm<int>> comm,
    Kokkos::View<DTK_CellTopology *, DeviceType>,
    Kokkos::View<unsigned int *, DeviceType>,
    Kokkos::View<double **, DeviceType>, Kokkos::View<double **, DeviceType>,
    Kokkos::View<LocalOrdinal *, DeviceType>,
    Kokkos::View<unsigned int *, DeviceType>,
    Kokkos::View<DTK_Quadrature *, DeviceType>, DTK_FEType )
    : _distributor( comm )
{
    DTK_INSIST( false );
}

template <typename DeviceType>
void Interpolation<DeviceType>::getReferencePoints(
    Kokkos::View<DataTransferKit::Point *, DeviceType> &phys_points,
    Kokkos::View<DataTransferKit::Point *, DeviceType> &ref_points,
    Kokkos::View<bool *, DeviceType> &pt_in_cell )
{
    // Merge block Kokkos::View in the same order as they were first
    // imported. This is necessary to reuse the distributor.
    unsigned int n_exported_points = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        n_exported_points += _reference_points[topo_id].extent( 0 );

    Kokkos::View<Point *, DeviceType> exported_ref_points(
        "exported_ref_points", n_exported_points );
    Kokkos::View<Point *, DeviceType> exported_phys_points(
        "exported_phys_points", n_exported_points );
    Kokkos::View<bool *, DeviceType> exported_pt_in_cell( "exported_pt_in_cell",
                                                          n_exported_points );
    unsigned int dim = _dim;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        Kokkos::View<unsigned int *, DeviceType> map =
            _point_indices_map[topo_id];
        Kokkos::View<double **, DeviceType> reference_points =
            _reference_points[topo_id];
        Kokkos::View<double **, DeviceType> phys_points =
            _filtered_points[topo_id];
        Kokkos::View<bool *, DeviceType> point_in_cell =
            _point_in_cell[topo_id];
        using ExecutionSpace = typename DeviceType::execution_space;
        unsigned int const n_topo_points = reference_points.extent( 0 );
        Kokkos::parallel_for(
            REGION_NAME( "merge_block_view" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_topo_points ),
            KOKKOS_LAMBDA( int const i ) {
                for ( unsigned int d = 0; d < dim; ++d )
                {
                    exported_ref_points( map( i ) )[d] =
                        reference_points( i, d );
                    exported_phys_points( map( i ) )[d] = phys_points( i, d );
                }
                exported_pt_in_cell( map( i ) ) = point_in_cell( i );
            } );
        Kokkos::fence();
    }

    // Communicate the results back to the calling processors
    unsigned int n_imports = _distributor.getTotalReceiveLength();
    Kokkos::realloc( phys_points, n_imports );
    Kokkos::realloc( ref_points, n_imports );
    Kokkos::realloc( pt_in_cell, n_imports );

    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _distributor, exported_phys_points, phys_points );
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _distributor, exported_ref_points, ref_points );
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _distributor, exported_pt_in_cell, pt_in_cell );
}

template <typename DeviceType>
void Interpolation<DeviceType>::convertMesh(
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<double **, DeviceType> coordinates )
{
    // First, we get the number of nodes for each topology to get initialize
    // _block_cells at the right size.
    // shards objects don't exist on the GPU, so we need to work on the CPU and
    // then copy the data back
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo(
        "n_nodes_per_topo" );
    auto n_nodes_per_topo_host = Kokkos::create_mirror_view( n_nodes_per_topo );
    for ( int i = 0; i < DTK_N_TOPO; ++i )
        n_nodes_per_topo_host( i ) = _cell_topologies[i].getNodeCount();
    Kokkos::deep_copy( n_nodes_per_topo, n_nodes_per_topo_host );

    // Create an array of Kokkos::View where each View is a block of cells with
    // the same topology
    for ( int i = 0; i < DTK_N_TOPO; ++i )
    {
        _block_cells[i] = Kokkos::View<double ***, DeviceType>(
            "block_cells_" + std::to_string( i ), _n_cells_per_topo[i],
            n_nodes_per_topo_host( i ), _dim );
    }

    // Compute the offset associated to each cell in the coordinates View.
    unsigned int const n_cells = _cell_topologies_view.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> nodes_per_cell( "nodes_per_cell",
                                                             n_cells );
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view =
        _cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> node_offset( "node_offset",
                                                          n_cells );
    computeNodeOffset( n_cells, cell_topo_view, n_nodes_per_topo,
                       nodes_per_cell, node_offset );

    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<unsigned int *, DeviceType> offset( "offset", n_cells );
    unsigned int const dim = _dim;
    Kokkos::View<DataTransferKit::Box *, DeviceType> bounding_boxes =
        _bounding_boxes;
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
                    computeBlockCellsBoundingBox(
                        dim, i, n_nodes_per_topo( topo_id ), node_offset( i ),
                        topo_id, cells, offset, coordinates, block_cells,
                        bounding_boxes, bounding_box_to_cell );
                }
            } );
        Kokkos::fence();
    }
}

template <typename DeviceType>
void Interpolation<DeviceType>::performDistributedSearch(
    Kokkos::View<double **, DeviceType> points_coord,
    Kokkos::View<DataTransferKit::Point *, DeviceType> &imported_points,
    Kokkos::View<int *, DeviceType> &imported_query_ids,
    Kokkos::View<int *, DeviceType> &imported_cell_indices )
{
    DataTransferKit::DistributedSearchTree<DeviceType> distributed_tree(
        _comm, _bounding_boxes );

    unsigned int const n_points = points_coord.extent( 0 );

    // Build the queries
    // FIXME do not use Within predicate because it requires a radius, i.e., it
    // is mesh dependent
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<DataTransferKit::Details::Within *, DeviceType> queries(
        "queries", n_points );
    Kokkos::parallel_for( "register_queries",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
                          KOKKOS_LAMBDA( int i ) {
                              queries( i ) = DataTransferKit::Details::within(
                                  {{points_coord( i, 0 ), points_coord( i, 1 ),
                                    points_coord( i, 2 )}},
                                  1e-14 );
                          } );
    Kokkos::fence();

    // Perform the distributed search
    Kokkos::View<int *, DeviceType> indices( "indices" );
    Kokkos::View<int *, DeviceType> offset( "offset" );
    Kokkos::View<int *, DeviceType> ranks( "ranks" );
    distributed_tree.query( queries, indices, offset, ranks );

    // Create the distributor and its reverse that will be used to communicate
    // the results back to the calling processors.
    auto ranks_host = Kokkos::create_mirror_view( ranks );
    Kokkos::deep_copy( ranks_host, ranks );
    Tpetra::Distributor distributed_search_distributor( _comm );
    unsigned int const n_imports =
        distributed_search_distributor.createFromSends(
            Teuchos::ArrayView<int const>( ranks_host.data(),
                                           ranks_host.extent( 0 ) ) );

    // Communicate cell indices
    Kokkos::realloc( imported_cell_indices, n_imports );
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        distributed_search_distributor, indices, imported_cell_indices );
    // Duplicate the points_coord for the communication. Duplicating the points
    // allows us to use the same distributor.
    unsigned int const indices_size = indices.extent( 0 );
    Kokkos::View<DataTransferKit::Point *, DeviceType> exported_points(
        "exported_points", indices_size );
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
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        distributed_search_distributor, exported_points, imported_points );

    // Communicate the query_ids. We communicate the query_ids to keep track of
    // which points is associated to which query.
    Kokkos::realloc( imported_query_ids, n_imports );
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        distributed_search_distributor, exported_query_ids,
        imported_query_ids );

    // We would like to get the reverse of the distributor but the following
    // code hangs during the communication
    //
    // Tpetra::Distributor
    // tmp_1(*(distributed_search_distributor.getReverse()));
    // Tpetra::Distributor tmp_2(*(tmp_1.getReverse()));
    // DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
    //     tmp_2, exported_points, imported_points );
    //
    // tmp_2 should be the same as distributed_search_distributor yet the code
    // hangs. So instead, create a new distributor from scratch. Note that I
    // cannot reproduce this simple on a simple test case.
    Kokkos::View<int *, DeviceType> own_rank( "own_rank", ranks.extent( 0 ) );
    Kokkos::deep_copy( own_rank, _comm->getRank() );
    Kokkos::View<int *, DeviceType> imported_ranks( "imported_ranks",
                                                    n_imports );
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        distributed_search_distributor, own_rank, imported_ranks );
    auto imported_ranks_host = Kokkos::create_mirror_view( imported_ranks );
    Kokkos::deep_copy( imported_ranks_host, imported_ranks );
    _distributor.createFromSends( Teuchos::ArrayView<int const>(
        imported_ranks_host.data(), imported_ranks_host.extent( 0 ) ) );
}

template <typename DeviceType>
void Interpolation<DeviceType>::filterTopology(
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    Kokkos::View<int *, DeviceType> cell_indices,
    Kokkos::View<DataTransferKit::Point *, DeviceType> points,
    Kokkos::View<int *, DeviceType> query_ids,
    Kokkos::View<unsigned int *, DeviceType> point_indices_map,
    Kokkos::View<int *, DeviceType> filtered_cell_indices,
    Kokkos::View<double **, DeviceType> filtered_points,
    Kokkos::View<int *, DeviceType> filtered_query_ids )
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
            }
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void Interpolation<DeviceType>::findReferencePoints(
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<double **, DeviceType> nodes_coordinates,
    Kokkos::View<double **, DeviceType> points_coordinates )
{
    // Initialize _bounding_box_to_cell to an invalid state
    Kokkos::deep_copy( _bounding_box_to_cell, static_cast<unsigned int>( -1 ) );

    // Compute the number of cells of each of the supported topologies.
    computeTopologies();

    // Convert the cells and nodes_coordinates View to _block_cells.
    convertMesh( cells, nodes_coordinates );

    // Perform the distributed search
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<DataTransferKit::Point *, DeviceType> imported_points(
        "imported_points", 0 );
    Kokkos::View<int *, DeviceType> imported_query_ids( "imported_query_ids",
                                                        0 );
    Kokkos::View<int *, DeviceType> imported_cell_indices( "imported_indices",
                                                           0 );
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
                                  imported_query_ids, imported_cell_indices );
    }
    else
        performDistributedSearch( points_coordinates, imported_points,
                                  imported_query_ids, imported_cell_indices );

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
    auto topo_size_host = Kokkos::create_mirror_view( topo_size );
    Kokkos::deep_copy( topo_size_host, topo_size );
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        _filtered_points[i] = Kokkos::View<double **, DeviceType>(
            "filtered_points_" + std::to_string( i ), topo_size_host( i ),
            _dim );
        _filtered_cell_indices[i] = Kokkos::View<int *, DeviceType>(
            "filtered_cell_indices_" + std::to_string( i ),
            topo_size_host( i ) );
        _filtered_query_ids[i] = Kokkos::View<int *, DeviceType>(
            "filtered_query_ids_" + std::to_string( i ), topo_size_host( i ) );
        _point_indices_map[i] = Kokkos::View<unsigned int *, DeviceType>(
            "point_indices_map_" + std::to_string( i ), topo_size_host( i ) );
        _reference_points[i] = Kokkos::View<double **, DeviceType>(
            "reference_points_" + std::to_string( i ), topo_size_host( i ),
            _dim );
        _point_in_cell[i] = Kokkos::View<bool *, DeviceType>(
            "point_in_cell_" + std::to_string( i ), topo_size_host( i ) );
    }

    // Check if the points are in the cells
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        if ( _block_cells[topo_id].extent( 0 ) != 0 )
            performPointInCell(
                _block_cells[topo_id], _cell_topologies[topo_id],
                imported_cell_indices, imported_points, imported_query_ids,
                topo, topo_id, _filtered_points[topo_id],
                _filtered_cell_indices[topo_id], _filtered_query_ids[topo_id],
                _point_indices_map[topo_id], _reference_points[topo_id],
                _point_in_cell[topo_id] );
}

template <typename DeviceType>
void Interpolation<DeviceType>::computeNodeOffset(
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
void Interpolation<DeviceType>::fillCellDofIds(
    unsigned int const topo_id, unsigned int const n_cells,
    unsigned int const n_nodes_per_topo,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view,
    Kokkos::View<unsigned int *, DeviceType> offset,
    Kokkos::View<unsigned int *, DeviceType> node_offset,
    Kokkos::View<LocalOrdinal *, DeviceType> object_dof_ids,
    Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for(
        REGION_NAME( "fill_cell_dof_ids" + std::to_string( topo_id ) ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
        KOKKOS_LAMBDA( int const i ) {
            if ( cell_topo_view( i ) == topo_id )
            {
                unsigned int const k = offset( i );
                for ( unsigned int node = 0; node < n_nodes_per_topo; ++node )
                {
                    unsigned int const n = node_offset( i ) + node;
                    cell_dofs_ids( k, node ) = object_dof_ids( n );
                }
            }
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void Interpolation<DeviceType>::computeTopologies()
{
    // This needs to be done on the host because we will use _n_cells_per_topo
    // to allocate Kokkos::View and this needs to be done on the host
    auto cell_topologies_host =
        Kokkos::create_mirror_view( _cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_host, _cell_topologies_view );
    unsigned int const n_cells = _cell_topologies_view.extent( 0 );
    DTK_CellTopology dtk_cell_topo;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        dtk_cell_topo = cell_topologies_host( i );
        _n_cells_per_topo[dtk_cell_topo] += 1;
    }

    // We do not support meshes that contains both 2D and 3D cells. All the
    // cells are either 2D or 3D
    _dim = _cell_topologies[dtk_cell_topo].getDimension();
#if HAVE_DTK_DBC
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        if ( _n_cells_per_topo[i] != 0 )
            DTK_REQUIRE( _cell_topologies[dtk_cell_topo].getDimension() ==
                         _dim );
    }
#endif

    // PointInCell does not support all Intrepid2 topologies
    DTK_REQUIRE( _n_cells_per_topo[DTK_TET_11] == 0 );
    DTK_REQUIRE( _n_cells_per_topo[DTK_HEX_20] == 0 );
    DTK_REQUIRE( _n_cells_per_topo[DTK_WEDGE_15] == 0 );
}

template <typename DeviceType>
void Interpolation<DeviceType>::performPointInCell(
    Kokkos::View<double ***, DeviceType> cells,
    shards::CellTopology const &cell_topology,
    Kokkos::View<int *, DeviceType> imported_cell_indices,
    Kokkos::View<DataTransferKit::Point *, DeviceType> imported_points,
    Kokkos::View<int *, DeviceType> imported_query_ids,
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    Kokkos::View<double **, DeviceType> filtered_points,
    Kokkos::View<int *, DeviceType> filtered_cell_indices,
    Kokkos::View<int *, DeviceType> filtered_query_ids,
    Kokkos::View<unsigned int *, DeviceType> point_indices_map,
    Kokkos::View<double **, DeviceType> reference_points,
    Kokkos::View<bool *, DeviceType> point_in_cell )
{
    filterTopology( topo, topo_id, imported_cell_indices, imported_points,
                    imported_query_ids, point_indices_map,
                    filtered_cell_indices, filtered_points,
                    filtered_query_ids );

    // Perform the PointInCell search
    DataTransferKit::PointInCell<DeviceType>::search(
        filtered_points, cells, filtered_cell_indices, cell_topology,
        reference_points, point_in_cell );
}
}

// Explicit instantiation macro
#define DTK_INTERPOLATION_INSTANT( NODE )                                      \
    template class Interpolation<typename NODE::device_type>;

#endif
