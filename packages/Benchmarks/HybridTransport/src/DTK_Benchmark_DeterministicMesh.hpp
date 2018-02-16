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
 * \file DTK_Benchmark_DeterministicMesh.hpp
 * \brief Deterministic mesh interface for the hybrid transport benchmark.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DETERMINISTICMESH_HPP
#define DTK_DETERMINISTICMESH_HPP

#include "CartesianMesh.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <vector>

namespace DataTransferKit
{
namespace Benchmark
{
//---------------------------------------------------------------------------//
/*!
 * \class DeterministicMesh
 * \brief Mesh and partitioner for deterministic transport simulations.
 *
 * This partitioner does a very simple partitioning where the client specifies
 * the number of partitions in \e (i,j).  The partitioner tries to make each
 * block (mesh on a processor) the same size.  If the number of cells in each
 * direction does not divide evenly, then cells are added to each direction
 * starting at \e (i=0,j=0). The mesh is decomposed into \e B blocks.
 *
 * The resulting partitioning is used to create the needed data in the base
 * class with mesh data being accessed through the base class interface.
 */
class DeterministicMesh : public CartesianMesh
{
  public:
    // Uniform cell size constructor. Will determine an "optimal" number of
    // blocks based on the communicator size.
    DeterministicMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                       const int num_cells_i, const int num_cells_j,
                       const int num_cells_k, const double delta_x,
                       const double delta_y, const double delta_z );

    // Uniform cell size/block size constructor.
    DeterministicMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                       const int num_i_blocks, const int num_j_blocks,
                       const int num_cells_i, const int num_cells_j,
                       const int num_cells_k, const double delta_x,
                       const double delta_y, const double delta_z );

    // Global/edge constructor.
    DeterministicMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                       const int num_i_blocks, const int num_j_blocks,
                       const std::vector<double> &global_x_edges,
                       const std::vector<double> &global_y_edges,
                       const std::vector<double> &global_z_edges );

  private:
    // Create uniform global edges and partition the mesh.
    void
    createEdgesAndPartition( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                             const int num_i_blocks, const int num_j_blocks,
                             const int num_cells_i, const int num_cells_j,
                             const int num_cells_k, const double delta_x,
                             const double delta_y, const double delta_z );

    // Partition the mesh.
    void partition( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                    const int num_i_blocks, const int num_j_blocks,
                    const std::vector<double> &global_x_edges,
                    const std::vector<double> &global_y_edges,
                    const std::vector<double> &global_z_edges );
};

//---------------------------------------------------------------------------//

} // end namespace Benchmark
} // end namespace DataTransferKit

#endif // end DTK_DETERMINISTICMESH_HPP
