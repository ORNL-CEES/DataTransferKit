/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_EXODUSASCIIGENERATOR_IMPL_HPP
#define DTK_EXODUSASCIIGENERATOR_IMPL_HPP

#include <Teuchos_Array.hpp>

#include <Tpetra_Distributor.hpp>

#include <iostream>
#include <fstream>
#include <limits>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template<class SourceDevice, class TargetDevice>
ExodusAsciiGenerator<SourceDevice,TargetDevice>::ExodusAsciiGenerator(
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const std::string& source_coord_file,
    const std::string& source_connectivity_file,
    const std::string& target_coord_file,
    const std::string& target_connectivity_file )
    : _comm( comm )
    , _src_coord_file( source_coord_file )
    , _src_connectivity_file( source_connectivity_file )
    , _tgt_coord_file( target_coord_file )
    , _tgt_connectivity_file( target_connectivity_file )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Create a problem where all points are uniquely owned (i.e. no ghosting)
template<class SourceDevice, class TargetDevice>
void ExodusAsciiGenerator<SourceDevice,TargetDevice>::createOneToOneProblem(
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,SourceDevice>& src_coords,
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,TargetDevice>& tgt_coords )
{
    // Partition the source coordinates in the x direction.
    partitionOneToOne( 0, _src_coord_file, src_coords );

    // Partition the target coordinates in the y direction.
    partitionOneToOne( 1, _tgt_coord_file, tgt_coords );
}

//---------------------------------------------------------------------------//
// Create a general problem where points max exist on multiple
// processors. Points have a unique global id.
template<class SourceDevice, class TargetDevice>
void ExodusAsciiGenerator<SourceDevice,TargetDevice>::createGeneralProblem(
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,SourceDevice>& src_coords,
    Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,SourceDevice>& src_gids,
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,TargetDevice>& tgt_coords,
    Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,TargetDevice>& tgt_gids )
{
    // Partition the source coordinates in the x direction.
    partitionGhostedConnectivity( 0, _src_coord_file, _src_connectivity_file,
                                  src_coords, src_gids );

    // Partition the target coordinates in the y direction.
    partitionGhostedConnectivity( 1, _tgt_coord_file, _tgt_connectivity_file,
                                  tgt_coords, tgt_gids );
}

//---------------------------------------------------------------------------//
template<class SourceDevice, class TargetDevice>
void ExodusAsciiGenerator<SourceDevice,TargetDevice>::getNodeDataFromFile(
    const std::string& coord_file,
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial>& host_coords,
    Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,Kokkos::Serial>& host_gids )
{
    // Only populate views on rank 0.
    if ( 0 == _comm->getRank() )
    {
        // Open the file.
        std::ifstream file;
        file.open( coord_file );

        // Get the number of nodes.
        int num_nodes = 0;
        file >> num_nodes;

        // Allocate the coordinate and global id arrray.
        host_coords = Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial>(
            "host_coords", num_nodes, 3 );
        host_gids = Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,Kokkos::Serial>(
            "host_gids", num_nodes );

        // Get the id and coordinate data.
        for ( int n = 0; n < num_nodes; ++n )
        {
            file >> host_gids(n);
            file >> host_coords(n,0);
            file >> host_coords(n,1);
            file >> host_coords(n,2);
        }

        // Close the file.
        file.close();
    }
    else
    {
        host_coords = Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial>(
            "host_coords", 0, 3 );
        host_gids = Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,Kokkos::Serial>(
            "host_gids", 0 );
    }
}

//---------------------------------------------------------------------------//
// Get the min and max of a given coordinate dimension.
template<class SourceDevice, class TargetDevice>
void ExodusAsciiGenerator<SourceDevice,TargetDevice>::dimensionMinAndMax(
    const int dim,
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial> coords,
    double& min,
    double& max )
{
    max = -std::numeric_limits<Coordinate>::max();
    min = std::numeric_limits<Coordinate>::max();
    int num_node = coords.extent(0);
    for ( int n = 0; n < num_node; ++n )
    {
        max = std::max( max, coords(n,dim) );
        min = std::min( min, coords(n,dim) );
    }
}

//---------------------------------------------------------------------------//
// Partition a point cloud in a given dimension with one-to-one mapping.
template<class SourceDevice, class TargetDevice>
template<class Device>
void ExodusAsciiGenerator<SourceDevice,TargetDevice>::partitionOneToOne(
    const int dim,
    const std::string& coord_file,
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Device>& device_coords )
{
    // Get the node data.
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial> host_coords;
    Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,Kokkos::Serial> host_gids;
    getNodeDataFromFile( coord_file, host_coords, host_gids );

    // Partition the coordinates in the right direction. Figure out the min
    // and max coordinates.
    Coordinate dim_max, dim_min;
    dimensionMinAndMax( dim, host_coords, dim_min, dim_max );

    // Build a communication plan.
    int num_node = host_coords.extent( 0 );
    Teuchos::Array<int> export_ranks( num_node );
    double dim_frac = 0.0;
    for ( int n = 0; n < num_node; ++n )
    {
        dim_frac = (host_coords(n,dim)-dim_min) / (dim_max-dim_min);
        export_ranks[n] = (1.0 == dim_frac)
                          ? _comm->getSize() - 1
                          : std::floor( dim_frac * _comm->getSize() );
    }
    Tpetra::Distributor distributor( _comm );
    int num_node_import = distributor.createFromSends( export_ranks() );

    // Send the sources to their new destination.
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial> import_coords(
        "import_coords", num_node_import, 3 );
    for ( int d = 0; d < 3; ++d )
    {
        distributor.doPostsAndWaits(
            Kokkos::subview(host_coords,Kokkos::ALL,d),
            1,
            Kokkos::subview(import_coords,Kokkos::ALL,d) );
    }

    // Move the coordinates to the device.
    device_coords =
        Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Device>( "device_coords", num_node_import, 3 );
    Kokkos::deep_copy( device_coords, import_coords );
}

//---------------------------------------------------------------------------//
// Partition a point cloud in a given dimension with ghosted connectivity
// mapping.
template<class SourceDevice, class TargetDevice>
template<class Device>
void ExodusAsciiGenerator<SourceDevice,TargetDevice>::partitionGhostedConnectivity(
    const int dim,
    const std::string& coord_file,
    const std::string& connectivity_file,
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Device>& device_coords,
    Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,Device> device_gids )
{
    // Get the node data.
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial> host_coords;
    Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,Kokkos::Serial> host_gids;
    getNodeDataFromFile( coord_file, host_coords, host_gids );

    // Partition the source coordinates in the dimension. Figure out the min
    // and max coordinates.
    Coordinate dim_max, dim_min;
    dimensionMinAndMax( dim, host_coords, dim_min, dim_max );

    // Load the connectivity file.
    std::ifstream conn_file;
    conn_file.open( connectivity_file );

    // Partition based on the dimension coordinate of the first node in each cell. All
    // nodes belonging to that cell will be sent to that rank.
    Teuchos::Array<GlobalOrdinal> export_gids(0);
    Teuchos::Array<int> export_ranks(0);
    Teuchos::Array<Coordinate> export_coords(0);

    // Read connectivity data on rank 0.
    if ( 0 == _comm->getRank() )
    {
        int num_cells = 0;
        conn_file >> num_cells;
        GlobalOrdinal cell_gid;
        GlobalOrdinal node_gid;
        int cell_num_node;
        std::string cell_topo;
        double x_frac = 0.0;
        int send_rank;
        for ( int c = 0; c < num_cells; ++c )
        {
            // Get the cell data.
            conn_file >> cell_gid;
            conn_file >> cell_topo;
            conn_file >> cell_num_node;

            // Get the first cell node.
            conn_file >> node_gid;

            // Partition.
            x_frac = (host_coords(node_gid-1,dim)-dim_min) / (dim_max-dim_min);
            send_rank = (1.0 == x_frac)
                        ? _comm->getSize() - 1
                        : std::floor( x_frac * _comm->getSize() );
            export_gids.push_back( node_gid );
            export_ranks.push_back( send_rank );
            export_coords.push_back( host_coords(node_gid-1,0) );
            export_coords.push_back( host_coords(node_gid-1,1) );
            export_coords.push_back( host_coords(node_gid-1,2) );

            // Add the rest of the cell nodes to that sending rank.
            for ( int n = 1; n < cell_num_node; ++n )
            {
                conn_file >> node_gid;
                export_gids.push_back( node_gid );
                export_ranks.push_back( send_rank );
                export_coords.push_back( host_coords(node_gid-1,0) );
                export_coords.push_back( host_coords(node_gid-1,1) );
                export_coords.push_back( host_coords(node_gid-1,2) );
            }
        }
    }

    // Build a communication plan for the sources.
    Tpetra::Distributor distributor( _comm );
    int num_import = distributor.createFromSends( export_ranks() );

    // Redistribute the sources.
    Teuchos::Array<GlobalOrdinal> import_gids( num_import );
    Teuchos::Array<Coordinate> import_coords( 3 * num_import );
    distributor.doPostsAndWaits( export_gids().getConst(), 1, import_gids() );
    distributor.doPostsAndWaits( export_coords().getConst(), 3, import_coords() );

    // Move the sources to the device.
    Kokkos::View<Coordinate*,Kokkos::Serial>
        host_import_gids( "host_import_gids", num_import );
    Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Kokkos::Serial>
        host_import_coords( "host_import_coords", num_import, 3 );
    for ( int n = 0; n < num_import; ++n )
    {
        host_import_gids(n) = import_gids[n];
        for ( int d = 0; d < 3; ++d )
            host_import_coords(n,d) = import_coords[ 3*n + d ];
    }
    device_coords =
        Kokkos::View<Coordinate**,Kokkos::LayoutLeft,Device>( "coords", num_import, 3 );
    device_gids =
        Kokkos::View<GlobalOrdinal*,Kokkos::LayoutLeft,Device>( "gids", num_import, 3 );
    Kokkos::deep_copy( device_coords, host_import_coords );
    Kokkos::deep_copy( device_gids, host_import_gids );

    // Close the source file.
    conn_file.close();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end  DTK_EXODUSASCIIGENERATOR_IMPL_HPP
