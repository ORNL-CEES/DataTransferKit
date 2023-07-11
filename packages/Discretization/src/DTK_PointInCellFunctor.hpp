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

#ifndef DTK_POINT_IN_CELL_FUNCTOR_HPP
#define DTK_POINT_IN_CELL_FUNCTOR_HPP

#include <Intrepid2_CellTools_Serial.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
namespace Functor
{
template <typename CellType, typename DeviceType>
class PointInCell
{
  public:
    PointInCell( double threshold,
                 Kokkos::View<double **, DeviceType> physical_points,
                 Kokkos::View<double ***, DeviceType> cells,
                 Kokkos::View<int *, DeviceType> coarse_search_output_cells,
                 Kokkos::View<double **, DeviceType> reference_points,
                 Kokkos::View<bool *, DeviceType> point_in_cell )
        : _threshold( threshold )
        , _physical_points( physical_points )
        , _cells( cells )
        , _coarse_search_output_cells( coarse_search_output_cells )
        , _reference_points( reference_points )
        , _point_in_cell( point_in_cell )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( unsigned int const i ) const
    {
        // Extract the indices computed by the coarse search
        int const cell_index = _coarse_search_output_cells( i );
        // Get the subviews corresponding the reference point (dim), the
        // physical point (dim), the current cell (nodes, dim)
        using ExecutionSpace = typename DeviceType::execution_space;
        Kokkos::View<double *, Kokkos::LayoutStride, ExecutionSpace> ref_point(
            _reference_points, i, Kokkos::ALL() );
        Kokkos::View<double *, Kokkos::LayoutStride, ExecutionSpace> phys_point(
            _physical_points, i, Kokkos::ALL() );
        Kokkos::View<double **, Kokkos::LayoutStride, ExecutionSpace> nodes(
            _cells, cell_index, Kokkos::ALL(), Kokkos::ALL() );

        // Compute the reference point and return true if the
        // point is inside the cell
        Intrepid2::Impl::CellTools::Serial::mapToReferenceFrame<
            typename CellType::basis_type>( ref_point, phys_point, nodes );
        _point_in_cell[i] =
            CellType::topo_type::checkPointInclusion( ref_point, _threshold );
    }

  private:
    double _threshold;
    Kokkos::View<double **, DeviceType> _physical_points;
    Kokkos::View<double ***, DeviceType> _cells;
    Kokkos::View<int *, DeviceType> _coarse_search_output_cells;
    Kokkos::View<double **, DeviceType> _reference_points;
    Kokkos::View<bool *, DeviceType> _point_in_cell;
};
} // namespace Functor
} // namespace DataTransferKit

#endif
