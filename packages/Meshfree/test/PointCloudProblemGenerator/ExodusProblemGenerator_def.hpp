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

#ifndef DTK_EXODUSPROBLEMGENERATOR_DEF_HPP
#define DTK_EXODUSPROBLEMGENERATOR_DEF_HPP

#include <ArborX.hpp>
#include <DTK_DBC.hpp>

#include <netcdf.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template <class Scalar, class SourceDevice, class TargetDevice>
ExodusProblemGenerator<Scalar, SourceDevice, TargetDevice>::
    ExodusProblemGenerator( MPI_Comm comm,
                            const std::string &source_exodus_file,
                            const std::string &target_exodus_file )
    : _comm( comm )
    , _src_exodus_file( source_exodus_file )
    , _tgt_exodus_file( target_exodus_file )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Create a problem where all points are uniquely owned (i.e. no ghosting)
template <class Scalar, class SourceDevice, class TargetDevice>
void ExodusProblemGenerator<Scalar, SourceDevice, TargetDevice>::
    createUniquelyOwnedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, SourceDevice> &src_field,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, TargetDevice> &tgt_field )
{
    // Partition the source coordinates in the x direction.
    partitionUniquelyOwned( 0, _src_exodus_file, src_coords );

    // Partition the target coordinates in the y direction.
    partitionUniquelyOwned( 1, _tgt_exodus_file, tgt_coords );

    // Allocate the fields and initialize to zero.
    auto num_src = src_coords.extent( 0 );
    src_field = Kokkos::View<Scalar **, Kokkos::LayoutLeft, SourceDevice>(
        "src_field", num_src, 1 );
    Kokkos::deep_copy( src_field, 0.0 );
    auto num_tgt = tgt_coords.extent( 0 );
    tgt_field = Kokkos::View<Scalar **, Kokkos::LayoutLeft, TargetDevice>(
        "tgt_field", num_tgt, 1 );
    Kokkos::deep_copy( tgt_field, 0.0 );
}

//---------------------------------------------------------------------------//
// Create a general problem where points may exist on multiple
// processors. Points have a unique global id.
template <class Scalar, class SourceDevice, class TargetDevice>
void ExodusProblemGenerator<Scalar, SourceDevice, TargetDevice>::
    createGhostedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, SourceDevice>
            &src_gids,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, SourceDevice> &src_field,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, TargetDevice>
            &tgt_gids,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, TargetDevice> &tgt_field )
{
    // Partition the source coordinates in the x direction.
    partitionGhostedConnectivity( 0, _src_exodus_file, src_coords, src_gids );

    // Partition the target coordinates in the y direction.
    partitionGhostedConnectivity( 1, _tgt_exodus_file, tgt_coords, tgt_gids );

    // Allocate the fields and initialize to zero.
    auto num_src = src_coords.extent( 0 );
    src_field = Kokkos::View<Scalar **, Kokkos::LayoutLeft, SourceDevice>(
        "src_field", num_src, 1 );
    Kokkos::deep_copy( src_field, 0.0 );
    auto num_tgt = tgt_coords.extent( 0 );
    tgt_field = Kokkos::View<Scalar **, Kokkos::LayoutLeft, TargetDevice>(
        "tgt_field", num_tgt, 1 );
    Kokkos::deep_copy( tgt_field, 0.0 );
}

//---------------------------------------------------------------------------//
// Read coordinate data from file and generate unqiue global ids for the
// points.
template <class Scalar, class SourceDevice, class TargetDevice>
template <class Device>
void ExodusProblemGenerator<Scalar, SourceDevice, TargetDevice>::
    getNodeDataFromFile(
        const std::string &exodus_file,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device> &coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device> &gids )
{
    // Only populate views on rank 0.
    int comm_rank;
    MPI_Comm_rank( _comm, &comm_rank );
    if ( comm_rank == 0 )
    {
        // Open the exodus file.
        int nc_id;
        DTK_CHECK_ERROR_CODE(
            nc_open( exodus_file.c_str(), NC_NOWRITE, &nc_id ) );

        // Get the number of nodes.
        auto num_nodes = getNetcdfDimensionLength( nc_id, "num_nodes" );

        // Allocate the coordinate and global id arrray.
        coords = Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device>(
            "coords", num_nodes, 3 );
        gids = Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device>(
            "gids", num_nodes );

        // Get the coordinate variable ids.
        int coord_var_id_x;
        DTK_CHECK_ERROR_CODE(
            nc_inq_varid( nc_id, "coordx", &coord_var_id_x ) );
        int coord_var_id_y;
        DTK_CHECK_ERROR_CODE(
            nc_inq_varid( nc_id, "coordy", &coord_var_id_y ) );
        int coord_var_id_z;
        DTK_CHECK_ERROR_CODE(
            nc_inq_varid( nc_id, "coordz", &coord_var_id_z ) );

        // Get the coordinates.
        DTK_CHECK_ERROR_CODE(
            nc_get_var_double( nc_id, coord_var_id_x, coords.data() ) );
        DTK_CHECK_ERROR_CODE( nc_get_var_double( nc_id, coord_var_id_y,
                                                 coords.data() + num_nodes ) );
        DTK_CHECK_ERROR_CODE( nc_get_var_double(
            nc_id, coord_var_id_z, coords.data() + 2 * num_nodes ) );

        // Close the exodus file.
        DTK_CHECK_ERROR_CODE( nc_close( nc_id ) );

        // Create unique global ids starting at 1.
        for ( size_t i = 0; i < num_nodes; ++i )
            gids( i ) = i + 1;
    }
    else
    {
        coords = Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device>(
            "coords", 0, 3 );
        gids = Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device>(
            "host_gids", 0 );
    }
}

//---------------------------------------------------------------------------//
// Partition a point cloud in a given dimension with one-to-one mapping.
template <class Scalar, class SourceDevice, class TargetDevice>
template <class Device>
void ExodusProblemGenerator<Scalar, SourceDevice, TargetDevice>::
    partitionUniquelyOwned( const int dim, const std::string &exodus_file,
                            Kokkos::View<Coordinate **, Kokkos::LayoutLeft,
                                         Device> &partitioned_coords )
{
    // Get the node data.
    Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device> export_coords;
    Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device> export_gids;
    getNodeDataFromFile( exodus_file, export_coords, export_gids );

    // Build a communication plan. Nodes are partitioned into equal spatial
    // bins in the given dimension. There is one spatial bin for each comm
    // rank.
    int num_node_export = export_coords.extent( 0 );
    Kokkos::View<int *, Kokkos::HostSpace> export_ranks( "export_proc_ids",
                                                         num_node_export );
    using ExecutionSpace = typename Device::execution_space;
    int comm_size;
    MPI_Comm_size( _comm, &comm_size );
    if ( 0 < num_node_export )
    {
        // Figure out the min and max coordinates in the given dimension.
        Coordinate dim_max, dim_min;
        std::tie( dim_min, dim_max ) = ArborX::minMax(
            ExecutionSpace{},
            Kokkos::subview( export_coords, Kokkos::ALL, dim ) );

        double dim_frac = 0.0;
        for ( int n = 0; n < num_node_export; ++n )
        {
            dim_frac =
                ( export_coords( n, dim ) - dim_min ) / ( dim_max - dim_min );
            export_ranks[n] = ( dim_frac < 1.0 )
                                  ? std::floor( dim_frac * comm_size )
                                  : comm_size - 1;
        }
    }
    ArborX::Details::Distributor<Device> distributor( _comm );
    int num_node_import =
        distributor.createFromSends( ExecutionSpace{}, export_ranks );

    // Send the coordinates to their new owning rank.
    partitioned_coords =
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device>(
            "partitioned_coords", num_node_import, 3 );
    ArborX::Details::DistributedTreeImpl<Device>::sendAcrossNetwork(
        ExecutionSpace{}, distributor, export_coords, partitioned_coords );
}

//---------------------------------------------------------------------------//
// Partition a point cloud in a given dimension with ghosted connectivity
// mapping.
template <class Scalar, class SourceDevice, class TargetDevice>
template <class Device>
void ExodusProblemGenerator<Scalar, SourceDevice, TargetDevice>::
    partitionGhostedConnectivity(
        const int dim, const std::string &exodus_file,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device>
            &partitioned_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device>
            &partitioned_gids )
{
    // Get the node data.
    Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device> input_coords;
    Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device> input_gids;
    getNodeDataFromFile( exodus_file, input_coords, input_gids );

    // Partition based on the dimension coordinate of the first node in each
    // cell. All nodes belonging to that cell will be sent to that rank to
    // simulate an element-based partitioning. This means nodes belonging to
    // elements on partition boundaries will exist on multiple ranks.
    std::vector<GlobalOrdinal> export_gids;
    std::vector<int> export_ranks;
    std::vector<Coordinate> export_coords;
    std::set<std::pair<int, GlobalOrdinal>> unique_exports;
    using ExecutionSpace = typename Device::execution_space;

    // Read connectivity data on rank 0.
    int comm_rank;
    MPI_Comm_rank( _comm, &comm_rank );
    int comm_size;
    MPI_Comm_size( _comm, &comm_size );
    if ( comm_rank == 0 )
    {
        // Figure out the min and max coordinates.
        Coordinate dim_max, dim_min;
        std::tie( dim_min, dim_max ) =
            ArborX::minMax( ExecutionSpace{},
                            Kokkos::subview( input_coords, Kokkos::ALL, dim ) );

        // Open the exodus file.
        int nc_id;
        DTK_CHECK_ERROR_CODE(
            nc_open( exodus_file.c_str(), NC_NOWRITE, &nc_id ) );

        // Get the number of element blocks.
        auto num_el_blks = getNetcdfDimensionLength( nc_id, "num_el_blk" );

        // Loop over blocks. Block ids start at 1.
        GlobalOrdinal node_gid;
        double x_frac = 0.0;
        int send_rank;
        for ( size_t b = 1; b < num_el_blks + 1; ++b )
        {
            // Get the number of elements in the block.
            std::string num_elem_dim_name =
                "num_el_in_blk" + std::to_string( b );
            auto num_elem =
                getNetcdfDimensionLength( nc_id, num_elem_dim_name );

            // Get the number of nodes per element in the block.
            std::string node_per_elem_dim_name =
                "num_nod_per_el" + std::to_string( b );
            auto node_per_elem =
                getNetcdfDimensionLength( nc_id, node_per_elem_dim_name );

            // Allocate a temporary view to load the block connectivity data.
            Kokkos::View<int **, Kokkos::LayoutLeft, Device> connectivity(
                "connectivity", node_per_elem, num_elem );

            // Get the connectivity.
            std::string conn_var_name = "connect" + std::to_string( b );
            int conn_var_id;
            DTK_CHECK_ERROR_CODE(
                nc_inq_varid( nc_id, conn_var_name.c_str(), &conn_var_id ) );
            DTK_CHECK_ERROR_CODE(
                nc_get_var_int( nc_id, conn_var_id, connectivity.data() ) );

            // Loop over elements in the block
            for ( size_t e = 0; e < num_elem; ++e )
            {
                // Get the first element node - connectivity indices start at 1.
                node_gid = connectivity( 0, e );

                // Partition in the given dimension. Each comm rank is
                // assigned a spatial bin along the given dimension. All
                // elements that have their first node in this spatial bin are
                // assigned to that comm rank.
                x_frac = ( input_coords( node_gid - 1, dim ) - dim_min ) /
                         ( dim_max - dim_min );
                send_rank = ( x_frac < 1.0 ) ? std::floor( x_frac * comm_size )
                                             : comm_size - 1;

                // Only add this node/rank combo if we haven't already. This
                // keeps us from sending the same node to the same rank more
                // than once.
                bool inserted = false;
                std::tie( std::ignore, inserted ) = unique_exports.insert(
                    std::make_pair( send_rank, node_gid ) );
                if ( inserted )
                {
                    export_gids.push_back( node_gid );
                    export_ranks.push_back( send_rank );
                    export_coords.push_back( input_coords( node_gid - 1, 0 ) );
                    export_coords.push_back( input_coords( node_gid - 1, 1 ) );
                    export_coords.push_back( input_coords( node_gid - 1, 2 ) );
                }

                // Add the rest of the cell nodes to that sending rank.
                for ( size_t n = 1; n < node_per_elem; ++n )
                {
                    // Get the next node id.
                    node_gid = connectivity( n, e );

                    // Only add this node/rank combo if we haven't
                    // already. This keeps us from sending the same node to
                    // the same rank more than once.
                    if ( !unique_exports.count(
                             std::make_pair( send_rank, node_gid ) ) )
                    {
                        export_gids.push_back( node_gid );
                        export_ranks.push_back( send_rank );
                        export_coords.push_back(
                            input_coords( node_gid - 1, 0 ) );
                        export_coords.push_back(
                            input_coords( node_gid - 1, 1 ) );
                        export_coords.push_back(
                            input_coords( node_gid - 1, 2 ) );
                        unique_exports.insert(
                            std::make_pair( send_rank, node_gid ) );
                    }
                }
            }
        }

        // Close the exodus file.
        DTK_CHECK_ERROR_CODE( nc_close( nc_id ) );
    }

    // Build a communication plan for the sources.
    ArborX::Details::Distributor<Device> distributor( _comm );
    using ExecutionSpace = typename Device::execution_space;
    int num_import = distributor.createFromSends(
        ExecutionSpace{}, Kokkos::View<int const *, Kokkos::HostSpace,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                              export_ranks.data(), export_ranks.size() ) );

    // Redistribute the sources.
    Kokkos::View<GlobalOrdinal *, Kokkos::HostSpace> import_gids( "",
                                                                  num_import );
    Kokkos::View<Coordinate * [3], Kokkos::HostSpace> import_coords(
        "", num_import );
    ArborX::Details::DistributedTreeImpl<Device>::sendAcrossNetwork(
        ExecutionSpace{}, distributor,
        Kokkos::View<GlobalOrdinal /*const*/ *, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            export_gids.data(), export_gids.size() ),
        import_gids );
    ArborX::Details::DistributedTreeImpl<Device>::sendAcrossNetwork(
        ExecutionSpace{}, distributor,
        Kokkos::View<Coordinate /*const*/ * [3], Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            export_coords.data(), export_coords.size() / 3 ),
        import_coords );

    // Move the sources to the device.
    partitioned_gids =
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device>(
            "import_gids", num_import );
    Kokkos::deep_copy( partitioned_gids, import_gids );
    partitioned_coords =
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device>(
            "import_coords", num_import, 3 );
    Kokkos::deep_copy( partitioned_coords, import_coords );
}

//---------------------------------------------------------------------------//
// Given a netcdf handle and a dimension name get the length of that
// dimension.
template <class Scalar, class SourceDevice, class TargetDevice>
size_t ExodusProblemGenerator<Scalar, SourceDevice, TargetDevice>::
    getNetcdfDimensionLength( const int nc_id, const std::string &dim_name )
{
    int dim_id;
    DTK_CHECK_ERROR_CODE( nc_inq_dimid( nc_id, dim_name.c_str(), &dim_id ) );
    size_t dim_len;
    DTK_CHECK_ERROR_CODE( nc_inq_dimlen( nc_id, dim_id, &dim_len ) );
    return dim_len;
}

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

#endif // end  DTK_EXODUSPROBLEMGENERATOR_DEF_HPP
