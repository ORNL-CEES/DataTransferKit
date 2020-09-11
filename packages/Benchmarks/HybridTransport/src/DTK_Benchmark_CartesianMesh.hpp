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
/*!
 * \file DTK_Benchmark_CartesianMesh.hpp
 * \brief Local Cartesian mesh interface for the hybrid transport benchmark.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CARTESIANMESH_HPP
#define DTK_CARTESIANMESH_HPP

#include "DTK_Types.h"

#include <Kokkos_Core.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <vector>

namespace DataTransferKit
{
namespace Benchmark
{
//---------------------------------------------------------------------------//
/*
 * \class CartesianMesh
 *
 * \brief A local description of a Cartesian mesh partition (or a block in
 * Exnihilo parlance).
 *
 * Both adjoint and deterministic transport algorithms in Exnihilo operate on
 * Cartesian grids of different decompositions and replication. This class
 * captures a common representation of these grids once they have been
 * partitioned.
 *
 * The grids are partitioned into sets and blocks. A set is one instance of
 * the entire global mesh and this instance can be subdivided amongst multiple
 * processors where the mesh partition owned by a single processor is a
 * block. Some examples:
 *
 * Example 1: A mesh partitioned into 1 block and 4 sets. This mesh will live
 * on a communicator of size 4. The mesh will be be fully replicated on each
 * of the processors in the communicator (i.e. it will have no domain
 * decomposition).
 *
 * Example 2: A mesh partitioned into 4 blocks and 1 set. This mesh will live
 * on a communicator of size 4. The mesh will be subdivided onto each of the
 * four processors in the communicator (i.e. it will have no replication).
 *
 * Example 3: A mesh partitioned into 2 blocks and 3 sets. This mesh will live
 * on a communicator of size 6. The mesh will be replicated 3 times and each
 * replica of the mesh will be subdivided onto 2 processors via domain
 * decomposition.
 *
 * There are two ways to access the description of the mesh through this
 * class:
 *
 * - The node coordinates and their global ids can be accessed - these nodes
 *   will have ghosting.
 *
 * - The cells can be accessed via their global ids and either a connectivity
 *   array or the coordinates of their centers - these cells will be uniquely
 *   owned.
 *
 */
class CartesianMesh
{
  public:
    /*
     * \brief Data constructor.
     *
     * \param comm The parallel communicator over which the mesh is built.
     *
     * \param set_id The set id of the local mesh.
     *
     * \param block_id The block id of the local mesh.
     *
     * \param x_global_num_node The global number of nodes in the x direction.
     *
     * \param y_global_num_node The global number of nodes in the y direction.
     *
     * \param z_global_num_node The global number of nodes in the z direction.
     *
     * \param x_edge_offset The starting index of the nodes on this process in
     * the x direction.
     *
     * \param y_edge_offset The starting index of the nodes on this process in
     * the y direction.
     *
     * \param z_edge_offset The starting index of the nodes on this process in
     * the z direction.
     *
     * \param local_x_edges The local edges of the mesh (node locations) in
     * the x direction.
     *
     * \param local_y_edges The local edges of the mesh (node locations) in
     * the y direction.
     *
     * \param local_z_edges The local edges of the mesh (node locations) in
     * the z direction.
     */
    CartesianMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                   const int set_id, const int block_id, const int num_i_blocks,
                   const int num_j_blocks, const int num_k_blocks,
                   const int x_global_num_node, const int y_global_num_node,
                   const int x_edge_offset, const int y_edge_offset,
                   const int z_edge_offset,
                   const std::vector<double> &local_x_edges,
                   const std::vector<double> &local_y_edges,
                   const std::vector<double> &local_z_edges );

    // Destructor.
    virtual ~CartesianMesh() = default;

    // Get the communicator for this mesh.
    Teuchos::RCP<const Teuchos::Comm<int>> comm() const { return _comm; }

    // Get the set id of this mesh.
    int setId() const { return _set_id; }

    // Get the block id of this mesh.
    int blockId() const { return _block_id; }

    // Number of sets.
    int numSets() const { return _comm->getSize() / numBlocks(); }

    // Number of blocks.
    int numBlocks() const
    {
        return _num_i_blocks * _num_j_blocks * _num_k_blocks;
    }

    // Number of I blocks.
    int numBlocksI() const { return _num_i_blocks; }

    // Number of J blocks.
    int numBlocksJ() const { return _num_j_blocks; }

    // Number of K blocks.
    int numBlocksK() const { return _num_k_blocks; }

    // Get the local node global ids.
    Kokkos::View<GlobalOrdinal *> localNodeGlobalIds() const
    {
        return _local_node_global_ids;
    }

    // Get the local node coordinates.
    Kokkos::View<Coordinate **> localNodeCoordinates() const
    {
        return _local_node_coords;
    }

    // Get the local cell ids.
    Kokkos::View<GlobalOrdinal *> localCellGlobalIds() const
    {
        return _local_cell_global_ids;
    }

    // Get the local cell connectivities.
    Kokkos::View<LocalOrdinal **> localCellConnectivity() const
    {
        return _local_cell_connectivity;
    }

    // Get the local cell center coordinates.
    Kokkos::View<Coordinate **> localCellCenterCoordinates() const
    {
        return _local_cell_center_coords;
    }

  private:
    // Communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> _comm;

    // Set id of this mesh.
    int _set_id;

    // Block id of this mesh.
    int _block_id;

    // Number of global I blocks.
    int _num_i_blocks;

    // Number of global J blocks.
    int _num_j_blocks;

    // Number of global K blocks.
    int _num_k_blocks;

    // Local node global ids.
    Kokkos::View<GlobalOrdinal *> _local_node_global_ids;

    // Local node coordinates.
    Kokkos::View<Coordinate **> _local_node_coords;

    // Local cell ids.
    Kokkos::View<GlobalOrdinal *> _local_cell_global_ids;

    // Get the local cell connectivities.
    Kokkos::View<LocalOrdinal **> _local_cell_connectivity;

    // Local cell center coordinates.
    Kokkos::View<Coordinate **> _local_cell_center_coords;
};

//---------------------------------------------------------------------------//

} // end namespace Benchmark
} // end namespace DataTransferKit

#endif // end DTK_CARTESIANMESH_HPP
