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

#ifndef DTK_SPLINE_OPERATOR_DECL_HPP
#define DTK_SPLINE_OPERATOR_DECL_HPP

#include <DTK_CompactlySupportedRadialBasisFunctions.hpp>
#include <DTK_MultivariatePolynomialBasis.hpp>
#include <DTK_PointCloudOperator.hpp>

#include <Tpetra_CrsMatrix.hpp>

#include <Thyra_LinearOpBase.hpp>

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <mpi.h>

namespace DataTransferKit
{

/**
 * This class implements a function reconstruction technique for arbitrary point
 * cloud based on a spline discretization. In this method, support
 * and subsequently the data transfer operator is constructed through solutions
 * to local least square kernels defined by compactly supported radial basis
 * functions.
 *
 * The class is templated on the DeviceType, the radial basis function
 * (Wendland<0>, Wendland<2>, Wendland<4>, Wendland<6>, Wu<2>, Wu<4>,
 * Buhmann<2>, Buhmann<3>, or Buhmann<4>) and polynonial basis (<Constant, DIM>,
 * <Linear, DIM>, or <Quadratic, DIM>).
 */
template <typename DeviceType,
          typename CompactlySupportedRadialBasisFunction = Wendland<0>,
          typename PolynomialBasis = MultivariatePolynomialBasis<Linear, 3>>
class SplineOperator : public PointCloudOperator<DeviceType>
{
    static_assert( std::is_same<PolynomialBasis,
                                MultivariatePolynomialBasis<Linear, 3>>::value,
                   "Only implemented for linear basis functions!" );
    using LO = int;
    using GO = long long;
    using NO = Kokkos::Compat::KokkosDeviceWrapperNode<
        typename DeviceType::execution_space>;
    using SC = Coordinate;

    using CrsMatrix = Tpetra::CrsMatrix<SC, LO, GO, NO>;
    using Map = Tpetra::Map<LO, GO, NO>;
    using Operator = Tpetra::Operator<SC, LO, GO, NO>;
    using Vector = Tpetra::MultiVector<SC, LO, GO, NO>;

  public:
    using device_type = DeviceType;
    using ExecutionSpace = typename DeviceType::execution_space;
    using polynomial_basis = PolynomialBasis;
    using radial_basis_function = CompactlySupportedRadialBasisFunction;

    SplineOperator(
        MPI_Comm comm,
        Kokkos::View<Coordinate const **, DeviceType> source_points,
        Kokkos::View<Coordinate const **, DeviceType> target_points );

    void
    apply( Kokkos::View<double const *, DeviceType> source_values,
           Kokkos::View<double *, DeviceType> target_values ) const override;

  private:
    MPI_Comm _comm;

    // Prolongation operator.
    Teuchos::RCP<const Operator> S;

    // Coefficient matrix polynomial component.
    Teuchos::RCP<const Operator> P;

    // Coefficient matrix basis component.
    Teuchos::RCP<const Operator> M;

    // Evaluation matrix polynomial component.
    Teuchos::RCP<const Operator> Q;

    // Evaluation matrix basis component.
    Teuchos::RCP<const Operator> N;

    // Coupling matrix
    Teuchos::RCP<const Thyra::LinearOpBase<SC>> _thyra_operator;

    // Source vector
    Teuchos::RCP<Vector> _source;

    // Destination vector
    Teuchos::RCP<Vector> _destination;

    // Source vector
    Teuchos::RCP<Thyra::MultiVectorBase<SC>> _thyra_X;

    // Destination vector
    Teuchos::RCP<Thyra::MultiVectorBase<SC>> _thyra_Y;

    Teuchos::RCP<Operator> buildPolynomialOperator(
        Teuchos::RCP<const Map> domain_map, Teuchos::RCP<const Map> range_map,
        Kokkos::View<Coordinate const **, DeviceType> points );

    Teuchos::RCP<Operator> buildBasisOperator(
        Teuchos::RCP<const Map> domain_map, Teuchos::RCP<const Map> range_map,
        Kokkos::View<Coordinate const **, DeviceType> source_points,
        Kokkos::View<Coordinate const **, DeviceType> target_points,
        int const knn );
};

} // end namespace DataTransferKit

#endif
