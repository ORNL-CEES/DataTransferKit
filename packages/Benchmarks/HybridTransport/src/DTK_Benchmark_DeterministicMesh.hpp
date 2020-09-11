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
 * \file DTK_Benchmark_DeterministicMesh.hpp
 * \brief Deterministic mesh interface for the hybrid transport benchmark.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DETERMINISTICMESH_HPP
#define DTK_DETERMINISTICMESH_HPP

#include "DTK_Benchmark_CartesianMesh.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <memory>
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
 * This partitioner does a very simple partitioning where the partitioner
 * calculates some "optimal" number of blocks which is most square.  The
 * partitioner tries to make each block (mesh on a processor) the same size.
 * If the number of cells in each direction does not divide evenly, then cells
 * are added to each direction starting at \e (i=0,j=0). The mesh is
 * decomposed into \e B blocks.
 *
 * The resulting partitioning is used to create the needed data in the base
 * class with mesh data being accessed through the base class
 * interface. Really think about this class as a decorator for the base class
 * which provides an initial partitioning.
 */
class DeterministicMesh
{
  public:
    /*!
     * \brief Uniform cell size constructor. The global list of node locations
     * will be composed from the given number of cells in uniform intervals.
     *
     * \param comm The parallel communicator over which the mesh is built.
     *
     * \param num_cells_i The global number of mesh cells in the X direction.
     *
     * \param num_cells_j The global number of mesh cells in the Y direction.
     *
     * \param num_cells_k The global number of mesh cells in the Z direction.
     *
     * \param delta_x The size of the mesh cells in the X direction.
     *
     * \param delta_y The size of the mesh cells in the Y direction.
     *
     * \param delta_z The size of the mesh cells in the Z direction.
     */
    DeterministicMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                       const int num_cells_i, const int num_cells_j,
                       const int num_cells_k, const double delta_x,
                       const double delta_y, const double delta_z );

    /*!
     * \brief Global edge constructor. A global list of node locations will be
     * used to create and partition the mesh.
     *
     * \param comm The parallel communicator over which the mesh is built.
     *
     * \param global_x_edges Global list of node locations in the X direction.
     *
     * \param global_y_edges Global list of node locations in the Y direction.
     *
     * \param global_z_edges Global list of node locations in the Z direction.
     */
    DeterministicMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                       const std::vector<double> &global_x_edges,
                       const std::vector<double> &global_y_edges,
                       const std::vector<double> &global_z_edges );

    /*!
     * \brief Get the local Cartesian mesh owned by this process.
     */
    std::shared_ptr<CartesianMesh> cartesianMesh() const
    {
        return _cartesian_mesh;
    }

  private:
    // Partition the mesh.
    void partition( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                    const std::vector<double> &global_x_edges,
                    const std::vector<double> &global_y_edges,
                    const std::vector<double> &global_z_edges );

  private:
    // The Cartesian mesh owned by this process.
    std::shared_ptr<CartesianMesh> _cartesian_mesh;
};

//---------------------------------------------------------------------------//

} // end namespace Benchmark
} // end namespace DataTransferKit

#endif // end DTK_DETERMINISTICMESH_HPP
