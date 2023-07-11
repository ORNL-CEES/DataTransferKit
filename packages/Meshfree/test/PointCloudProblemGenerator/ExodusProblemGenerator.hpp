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

#ifndef DTK_EXODUSPROBLEMGENERATOR_HPP
#define DTK_EXODUSPROBLEMGENERATOR_HPP

#include "DTK_ConfigDefs.hpp"
#include "DTK_Types.h"
#include "PointCloudProblemGenerator.hpp"

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <functional>
#include <string>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Generate point cloud problem by reading exodus files.
//
// The generator reads exodus files and extracts the node coordinates and
// partitions them across the given communicator.
//
// Source files are partitioned in x where each rank gets an even subdivision
// of space in the x dimension:
//
// |          | o        |        |                      |
// |  o     o |       o  |  o  o  |  o o       o         |
// |    o     |  o     o |      o |        o      o      |
// |   o   o  |     o    | o      |     o    o  o        |
// |          |          |    o   |                      |
// | rank = 0 | rank = 1 | ...... | rank = comm_size - 1 |
//
// Target files are partitioned in y where each rank gets an even subdivision
// of space in the y dimension.
//
// ---------------------------------
//             o o   o     o
// rank = 0       o     o   o
//              o   o   o o
// ---------------------------------
//               o   o    o   o
// rank = 1    o  o   o o    o
//              o    o oo    o o
// ---------------------------------
//               oo  o o o o  oo
// ....       o   o  o  o  o  o
//              o  oo  o o oo o o
// ---------------------------------
//                   o o o oo
// rank =             o   o o o
//   comm_size - 1   o   oo  ooooo
//
// ---------------------------------
//
// In the uniquely owned case, mesh nodes are given one single, unique
// destination in the partitioning.
//
// In the ghosted case, mesh elements are given one unique destination and
// all nodes belonging to that element are sent to that
// destination. Therefore, nodes owned by multiple elements with different
// destinations will be sent to multiple ranks and thus ghosted.
//
template <class Scalar, class SourceDevice, class TargetDevice>
class ExodusProblemGenerator
    : public PointCloudProblemGenerator<Scalar, SourceDevice, TargetDevice>
{
  public:
    // Constructor.
    ExodusProblemGenerator( MPI_Comm comm,
                            const std::string &source_exodus_file,
                            const std::string &target_exodus_file );

    // Create a problem where all points are uniquely owned (i.e. no
    // ghosting). Both source and target fields have one component and are
    // initialized to zero.
    void createUniquelyOwnedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, SourceDevice> &src_field,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, TargetDevice> &tgt_field )
        override;

    // Create a general problem where points may exist on multiple
    // processors. Both source and target fields have 1 component and are
    // initialized to zero.
    void createGhostedProblem(
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
        override;

  private:
    // Get host views of node data from file.
    template <class Device>
    void getNodeDataFromFile(
        const std::string &exodus_file,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device> &coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device> &gids );

    // Partition a point cloud in a given dimension with one-to-one mapping.
    template <class Device>
    void partitionUniquelyOwned( const int dim, const std::string &exodus_file,
                                 Kokkos::View<Coordinate **, Kokkos::LayoutLeft,
                                              Device> &partitioned_coords );

    // Partition a point cloud in a given dimension with ghosted connectivity
    // mapping.
    template <class Device>
    void partitionGhostedConnectivity(
        const int dim, const std::string &exodus_file,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device>
            &partitioned_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device>
            &partitioned_gids );

    // Given a netcdf handle and a dimension name get the length of that
    // dimension.
    size_t getNetcdfDimensionLength( const int nc_id,
                                     const std::string &dim_name );

  private:
    // Comm
    MPI_Comm _comm;

    // Filenames
    std::string _src_exodus_file;
    std::string _tgt_exodus_file;
};

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes
//---------------------------------------------------------------------------//

#include "ExodusProblemGenerator_def.hpp"

//---------------------------------------------------------------------------//

#endif // end  DTK_EXODUSPROBLEMGENERATOR_HPP
