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

#ifndef DTK_MOVING_LEAST_SQUARES_OPERATOR_DECL_HPP
#define DTK_MOVING_LEAST_SQUARES_OPERATOR_DECL_HPP

#include <DTK_CompactlySupportedRadialBasisFunctions.hpp>
#include <DTK_MultivariatePolynomialBasis.hpp>
#include <DTK_PointCloudOperator.hpp>

#include <mpi.h>

namespace DataTransferKit
{

template <typename DeviceType,
          typename CompactlySupportedRadialBasisFunction = Wendland<0>,
          typename PolynomialBasis = MultivariatePolynomialBasis<Linear, 3>>
class MovingLeastSquaresOperator : public PointCloudOperator<DeviceType>
{
    using ExecutionSpace = typename DeviceType::execution_space;

  public:
    MovingLeastSquaresOperator(
        MPI_Comm comm,
        Kokkos::View<Coordinate const **, DeviceType> source_points,
        Kokkos::View<Coordinate const **, DeviceType> target_points );

    void
    apply( Kokkos::View<double const *, DeviceType> source_values,
           Kokkos::View<double *, DeviceType> target_values ) const override;

  private:
    MPI_Comm _comm;
    unsigned int const _n_source_points;
    Kokkos::View<int *, DeviceType> _offset;
    Kokkos::View<int *, DeviceType> _ranks;
    Kokkos::View<int *, DeviceType> _indices;
    Kokkos::View<double *, DeviceType> _coeffs;
};

} // end namespace DataTransferKit

#endif
