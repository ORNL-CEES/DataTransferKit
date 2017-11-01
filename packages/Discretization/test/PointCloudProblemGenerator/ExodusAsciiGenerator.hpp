/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_EXODUSASCIIGENERATOR_HPP
#define DTK_EXODUSASCIIGENERATOR_HPP

#include "DTK_ConfigDefs.hpp"
#include "DTK_Types.h"
#include "PointCloudProblemGenerator.hpp"

#include <Kokkos_View.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <string>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template <class SourceDevice, class TargetDevice>
class ExodusAsciiGenerator
    : public PointCloudProblemGenerator<SourceDevice, TargetDevice>
{
  public:
    // Constructor.
    ExodusAsciiGenerator( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                          const std::string &source_coord_file,
                          const std::string &source_connectivity_file,
                          const std::string &target_coord_file,
                          const std::string &target_connectivity_file );

    // Create a problem where all points are uniquely owned (i.e. no ghosting)
    void createUniquelyOwnedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords ) override;

    // Create a general problem where points may exist on multiple
    // processors. Points have a unique global id.
    void createGhostedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, SourceDevice>
            &src_gids,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, TargetDevice>
            &tgt_gids ) override;

  private:
    // Get host views of node data from file.
    void getNodeDataFromFile( const std::string &coord_file,
                              Kokkos::View<Coordinate **, Kokkos::LayoutLeft,
                                           Kokkos::Serial> &host_coords,
                              Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft,
                                           Kokkos::Serial> &host_gids );

    // Get the min and max of a given coordinate dimension.
    void dimensionMinAndMax(
        const int dim,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Kokkos::Serial> coords,
        double &min, double &max );

    // Partition a point cloud in a given dimension with one-to-one mapping.
    template <class Device>
    void partitionUniquelyOwned( const int dim, const std::string &coord_file,
                                 Kokkos::View<Coordinate **, Kokkos::LayoutLeft,
                                              Device> &device_coords );

    // Partition a point cloud in a given dimension with ghosted connectivity
    // mapping.
    template <class Device>
    void partitionGhostedConnectivity(
        const int dim, const std::string &coord_file,
        const std::string &connectivity_file,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, Device> &device_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, Device>
            &device_gids );

  private:
    // Comm
    Teuchos::RCP<const Teuchos::Comm<int>> _comm;

    // Filenames
    std::string _src_coord_file;
    std::string _src_connectivity_file;
    std::string _tgt_coord_file;
    std::string _tgt_connectivity_file;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes
//---------------------------------------------------------------------------//

#include "ExodusAsciiGenerator_def.hpp"

//---------------------------------------------------------------------------//

#endif // end  DTK_EXODUSASCIIGENERATOR_HPP
