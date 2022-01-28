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

#ifndef DTK_INTERPOLATION_DEF_HPP
#define DTK_INTERPOLATION_DEF_HPP

#include <DTK_FE.hpp>
#include <DTK_PointInCell.hpp>

namespace DataTransferKit
{
template <typename DeviceType>
Interpolation<DeviceType>::Interpolation(
    MPI_Comm comm, Mesh<DeviceType> const &mesh,
    Kokkos::View<Coordinate **, DeviceType> points_coordinates,
    Kokkos::View<LocalOrdinal *, DeviceType> cell_dof_ids, DTK_FEType fe_type )
    : _point_search( comm, mesh, points_coordinates )
{
    // Fill up _finite_element, i.e., fill up a map between topo_id and FE
    Topologies topologies;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        _finite_elements[topo_id] = getFE( topologies[topo_id].topo, fe_type );

    // Change the format of cell_dofs_ids
    filter_dofs_ids( mesh.cell_topologies, cell_dof_ids, fe_type );
}

template <typename DeviceType>
void Interpolation<DeviceType>::filter_dofs_ids(
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<LocalOrdinal *, DeviceType> cell_dof_ids, DTK_FEType fe_type )
{
    // We need to filter the dof_ids and only keep the cells where a point
    // was found. Because multiple points may be in the same cells, the
    // cells may be duplicated.
    // TODO do this on the device
    auto cell_topologies_host = Kokkos::create_mirror_view( cell_topologies );
    Kokkos::deep_copy( cell_topologies_host, cell_topologies );

    unsigned int const n_cells = cell_topologies.extent( 0 );
    // We need to compute the number of basis function for each cell because the
    // number of basis functions is different for HGRAD, HDIV, and HCURL.
    // Therefore, knowing the number of nodes in the topology is not enough.
    std::vector<unsigned int> dof_offset( n_cells + 1 );
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        auto fe = getFE( cell_topologies_host( i ), fe_type );
        dof_offset[i + 1] = dof_offset[i] + getCardinality<DeviceType>( fe );
    }

    auto cell_dof_ids_host = Kokkos::create_mirror_view( cell_dof_ids );
    Kokkos::deep_copy( cell_dof_ids_host, cell_dof_ids );
    std::array<std::vector<std::vector<unsigned int>>, DTK_N_TOPO>
        filtered_dof_ids;
    // For each topo_id (finite element type) we reformat cell_dof_ids
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int const n_dofs_per_cell =
            getCardinality<DeviceType>( _finite_elements[topo_id] );
        auto cell_indices_host =
            Kokkos::create_mirror_view( _point_search._cell_indices[topo_id] );
        Kokkos::deep_copy( cell_indices_host,
                           _point_search._cell_indices[topo_id] );

        // For each cell which contains a target point, we reformat cell_dof_ids
        for ( unsigned int i = 0;
              i < _point_search._query_ids[topo_id].extent( 0 ); ++i )
        {
            unsigned int const cell_id =
                _point_search
                    ._cell_indices_map[topo_id][cell_indices_host( i )];
            unsigned int const offset = dof_offset[cell_id];
            std::vector<unsigned int> current_cell_dof_ids( n_dofs_per_cell );
            for ( unsigned int j = 0; j < n_dofs_per_cell; ++j )
                current_cell_dof_ids[j] = cell_dof_ids_host( offset + j );
            filtered_dof_ids[topo_id].push_back( current_cell_dof_ids );
        }
    }

    // Copy in a Kokkos::View and then move it to the device
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int const fe_n_cells = filtered_dof_ids[topo_id].size();
        unsigned int const n_dofs_per_cell =
            ( fe_n_cells > 0 ) ? filtered_dof_ids[topo_id][0].size() : 0;
        _dofs_ids[topo_id] = Kokkos::View<LocalOrdinal **, DeviceType>(
            "cell_dofs_ids_" + std::to_string( topo_id ), fe_n_cells,
            n_dofs_per_cell );
        auto dofs_ids_host = Kokkos::create_mirror_view( _dofs_ids[topo_id] );
        for ( unsigned int i = 0; i < fe_n_cells; ++i )
            for ( unsigned int j = 0; j < n_dofs_per_cell; ++j )
                dofs_ids_host( i, j ) = filtered_dof_ids[topo_id][i][j];
        Kokkos::deep_copy( _dofs_ids[topo_id], dofs_ids_host );
    }
}

} // namespace DataTransferKit

// Explicit instantiation macro
#define DTK_INTERPOLATION_INSTANT( NODE )                                      \
    template class Interpolation<typename NODE::device_type>;

#endif
