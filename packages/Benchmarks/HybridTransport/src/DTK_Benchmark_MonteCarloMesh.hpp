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
 * \file DTK_Benchmark_MonteCarloMesh.hpp
 * \brief MonteCarlo mesh interface for the hybrid transport benchmark.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MONTECARLOMESH_HPP
#define DTK_MONTECARLOMESH_HPP

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
 * \class MonteCarloMesh
 * \brief Mesh and partitioner for Monte Carlo transport simulations.
 *
 * Monte Carlo grids are typically subdivided into a small number of "blocks",
 * which can be interpreted as grid partitions, and "sets" which can be
 * interpreted as a number of replications of the grid.
 *
 * The block decomposition, or spatial partitioning of the mesh, is always
 * computed with user provided planes indicating the spatial boundaries of
 * each block. The underlying grid is then subdivided into blocks based on the
 * plane locations. If a plane intersects a cell, that cell is assigned to
 * both blocks adjacent to the plane.
 *
 * Note that the boundary mesh planes don't have to contain the grid - we will
 * just use them to carve out partitions from the background grid.
 *
 * The resulting partitioning is used to create the needed data in the base
 * class with mesh data being accessed through the base class
 * interface. Really think about this class as a decorator for the base class
 * which provides an initial partitioning.
 */
class MonteCarloMesh
{
  public:
    /*!
     * \brief Uniform cell size constructor.
     *
     * Builds a regular grid with a uniform cell size and partitions it using
     * the input boundary mesh arrays.
     *
     * \param comm The communicator over which the mesh will be partitioned.
     *
     * \param num_sets The number of sets over which to decompose the
     * mesh. This is the total number of times the mesh is replicated.
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
     *
     * \param x_bnd_mesh The boundary mesh node locations in the x direction.
     *
     * \param y_bnd_mesh The boundary mesh node locations in the y direction.
     *
     * \param z_bnd_mesh The boundary mesh node locations in the z direction.
     */
    MonteCarloMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                    const int num_sets, const int num_cells_i,
                    const int num_cells_j, const int num_cells_k,
                    const double delta_x, const double delta_y,
                    const double delta_z, const std::vector<double> &x_bnd_mesh,
                    const std::vector<double> &y_bnd_mesh,
                    const std::vector<double> &z_bnd_mesh );

    /*!
     * \brief Global edge constructor. A global list of node locations will be
     * used to create the mesh. The mesh will be partitioned using the input
     * boundary mesh arrays.
     *
     * \param comm The parallel communicator over which the mesh is built.
     *
     * \param num_sets The number of sets over which to decompose the
     * mesh. This is the total number of times the mesh is replicated.
     *
     * \param global_x_edges Global list of node locations in the X direction.
     *
     * \param global_y_edges Global list of node locations in the Y direction.
     *
     * \param global_z_edges Global list of node locations in the Z direction.
     *
     * \param x_bnd_mesh The boundary mesh node locations in the x direction.
     *
     * \param y_bnd_mesh The boundary mesh node locations in the y direction.
     *
     * \param z_bnd_mesh The boundary mesh node locations in the z direction.
     */
    MonteCarloMesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                    const int num_sets,
                    const std::vector<double> &global_x_edges,
                    const std::vector<double> &global_y_edges,
                    const std::vector<double> &global_z_edges,
                    const std::vector<double> &x_bnd_mesh,
                    const std::vector<double> &y_bnd_mesh,
                    const std::vector<double> &z_bnd_mesh );

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
                    const int num_sets,
                    const std::vector<double> &global_x_edges,
                    const std::vector<double> &global_y_edges,
                    const std::vector<double> &global_z_edges,
                    const std::vector<double> &x_bnd_mesh,
                    const std::vector<double> &y_bnd_mesh,
                    const std::vector<double> &z_bnd_mesh );

    // Calculate local edge arrays.
    void computeLocalEdges( const std::vector<double> &global_edges,
                            const std::vector<double> &bnd_mesh,
                            const int my_block,
                            std::vector<double> &local_edges,
                            int &offset ) const;

  private:
    // The Cartesian mesh owned by this process.
    std::shared_ptr<CartesianMesh> _cartesian_mesh;
};

//---------------------------------------------------------------------------//

} // end namespace Benchmark
} // end namespace DataTransferKit

#endif // end DTK_MONTECARLOMESH_HPP
