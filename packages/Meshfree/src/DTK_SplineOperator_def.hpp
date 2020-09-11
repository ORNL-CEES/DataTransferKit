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

#ifndef DTK_SPLINE_OPERATOR_DEF_HPP
#define DTK_SPLINE_OPERATOR_DEF_HPP

#include <ArborX.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <DTK_DBC.hpp>
#include <DTK_DetailsMovingLeastSquaresOperatorImpl.hpp>
#include <DTK_DetailsNearestNeighborOperatorImpl.hpp> // fetch
#include <DTK_DetailsPolynomialMatrix.hpp>
#include <DTK_DetailsSplineProlongationOperator.hpp>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Thyra_DefaultAddedLinearOp.hpp>
#include <Thyra_DefaultMultipliedLinearOp.hpp>
#include <Thyra_DefaultScaledAdjointLinearOp.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

namespace DataTransferKit
{

template <typename DeviceType, typename CompactlySupportedRadialBasisFunction,
          typename PolynomialBasis>
Teuchos::RCP<
    typename SplineOperator<DeviceType, CompactlySupportedRadialBasisFunction,
                            PolynomialBasis>::Operator>
SplineOperator<DeviceType, CompactlySupportedRadialBasisFunction,
               PolynomialBasis>::
    buildBasisOperator(
        Teuchos::RCP<const Map> domain_map, Teuchos::RCP<const Map> range_map,
        Kokkos::View<Coordinate const **, DeviceType> source_points,
        Kokkos::View<Coordinate const **, DeviceType> target_points,
        int const knn )
{
    auto teuchos_comm = domain_map->getComm();
    auto teuchos_mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( teuchos_comm );
    MPI_Comm comm = ( *teuchos_mpi_comm->getRawMpiComm() )();

    int const num_source_points = source_points.extent( 0 );
    int const num_points = target_points.extent( 0 );

    ArborX::DistributedSearchTree<DeviceType> distributed_tree( comm,
                                                                source_points );
    DTK_CHECK( !distributed_tree.empty() );

    // Perform the actual search.
    auto queries =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::makeKNNQueries(
            target_points, knn );

    Kokkos::View<int *, DeviceType> offset( "offset", 0 );
    Kokkos::View<int *, DeviceType> ranks( "ranks", 0 );
    Kokkos::View<int *, DeviceType> indices( "indices", 0 );
    distributed_tree.query( queries, indices, offset, ranks );

    // Retrieve the coordinates of all points that met the predicates.
    auto source_points_with_halo =
        Details::NearestNeighborOperatorImpl<DeviceType>::fetch(
            comm, ranks, indices, source_points );

    auto transformed_source_points = Details::MovingLeastSquaresOperatorImpl<
        DeviceType>::transformSourceCoordinates( source_points_with_halo,
                                                 offset, target_points );

    // To build the radial basis function, we need to define the radius of
    // the radial basis function. Since we use kNN, we need to compute the
    // radius. We only need the coordinates of the source points because of
    // the transformation of the coordinates.
    auto radius =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::computeRadius(
            transformed_source_points, offset );

    // Build phi (weight matrix)
    auto phi =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::computeWeights(
            transformed_source_points, radius,
            CompactlySupportedRadialBasisFunction() );

    // Compute some helper arrays
    int comm_size = teuchos_comm->getSize();

    std::vector<int> cumulative_points_per_process( comm_size + 1 );
    MPI_Allgather( &num_source_points, 1, MPI_INT,
                   cumulative_points_per_process.data(), 1, MPI_INT, comm );

    // Exclusive prefix scan
    int total = 0;
    for ( int i = 0; i < comm_size; ++i )
    {
        int current_value = cumulative_points_per_process[i];
        cumulative_points_per_process[i] = total;
        total += current_value;
    }

    // Build matrix
    auto row_map = range_map;
    auto crs_matrix = Teuchos::rcp( new CrsMatrix( row_map, knn ) );

    for ( LO i = 0; i < num_points; ++i )
        for ( int j = offset( i ); j < offset( i + 1 ); ++j )
        {
            const auto row_index = row_map->getGlobalElement( i );
            const auto col_index =
                cumulative_points_per_process[ranks( j )] + indices( j );
            const auto value = phi( j );

            crs_matrix->insertGlobalValues( row_index,
                                            Teuchos::tuple<GO>( col_index ),
                                            Teuchos::tuple<SC>( value ) );
        }

    crs_matrix->fillComplete( domain_map, range_map );
    DTK_ENSURE( crs_matrix->isFillComplete() );

    return crs_matrix;
}

template <typename DeviceType, typename CompactlySupportedRadialBasisFunction,
          typename PolynomialBasis>
Teuchos::RCP<
    typename SplineOperator<DeviceType, CompactlySupportedRadialBasisFunction,
                            PolynomialBasis>::Operator>
SplineOperator<DeviceType, CompactlySupportedRadialBasisFunction,
               PolynomialBasis>::
    buildPolynomialOperator(
        Teuchos::RCP<const Map> domain_map, Teuchos::RCP<const Map> range_map,
        Kokkos::View<Coordinate const **, DeviceType> points )
{
    const int spatial_dim = points.extent( 1 );

    DTK_REQUIRE( spatial_dim == 3 );

    auto v = Details::MovingLeastSquaresOperatorImpl<
        DeviceType>::computeVandermonde2( points, PolynomialBasis() );

    return Teuchos::rcp(
        new PolynomialMatrix<SC, LO, GO, NO>( v, domain_map, range_map ) );
}

template <typename DeviceType, typename CompactlySupportedRadialBasisFunction,
          typename PolynomialBasis>
SplineOperator<DeviceType, CompactlySupportedRadialBasisFunction,
               PolynomialBasis>::
    SplineOperator(
        MPI_Comm comm,
        Kokkos::View<Coordinate const **, DeviceType> source_points,
        Kokkos::View<Coordinate const **, DeviceType> target_points )
    : _comm( comm )
{
    DTK_REQUIRE( source_points.extent_int( 1 ) ==
                 target_points.extent_int( 1 ) );
    // FIXME for now let's assume 3D
    DTK_REQUIRE( source_points.extent_int( 1 ) == 3 );
    constexpr int spatial_dim = 3;

    constexpr int knn = PolynomialBasis::size;

    // Step 0: build source and target maps
    auto teuchos_comm = Teuchos::rcp( new Teuchos::MpiComm<int>( comm ) );
    auto source_map = Teuchos::rcp(
        new Map( Teuchos::OrdinalTraits<GO>::invalid(),
                 source_points.extent( 0 ), 0 /*indexBase*/, teuchos_comm ) );
    auto target_map = Teuchos::rcp(
        new Map( Teuchos::OrdinalTraits<GO>::invalid(),
                 target_points.extent( 0 ), 0 /*indexBase*/, teuchos_comm ) );

    // Step 1: build matrices
    GO prolongation_offset = teuchos_comm->getRank() ? 0 : spatial_dim + 1;
    S = Teuchos::rcp( new SplineProlongationOperator<SC, LO, GO, NO>(
        prolongation_offset, source_map ) );
    auto prolongation_map = S->getRangeMap();

    // Build distributed search tree over the source points.
    // NOTE: M is not the M from the paper, but an extended size block
    // matrix
    M = buildBasisOperator( prolongation_map, prolongation_map, source_points,
                            source_points, knn );
    P = buildPolynomialOperator( prolongation_map, prolongation_map,
                                 source_points );
    N = buildBasisOperator( prolongation_map, target_map, source_points,
                            target_points, knn );
    Q = buildPolynomialOperator( prolongation_map, target_map, target_points );

    // Step 3: build Thyra operator: A = (Q + N)*[(P + M + P^T)^-1]*S
    auto thyraWrapper = []( Teuchos::RCP<const Operator> &op ) {
        auto thyra_range_vector_space =
            Thyra::createVectorSpace<SC>( op->getRangeMap() );
        auto thyra_domain_vector_space =
            Thyra::createVectorSpace<SC>( op->getDomainMap() );
        using ThyraOperator = Thyra::TpetraLinearOp<SC, LO, GO, NO>;
        auto thyra_op = Teuchos::rcp( new ThyraOperator() );
        Teuchos::rcp_const_cast<ThyraOperator>( thyra_op )
            ->constInitialize( thyra_range_vector_space,
                               thyra_domain_vector_space, op );
        return thyra_op;
    };

    auto thyra_S = thyraWrapper( S );
    auto thyra_M = thyraWrapper( M );
    auto thyra_N = thyraWrapper( N );
    auto thyra_P = thyraWrapper( P );
    auto thyra_Q = thyraWrapper( Q );

    // Create a transpose of P.
    Teuchos::RCP<const Thyra::LinearOpBase<SC>> thyra_P_T =
        Thyra::transpose<SC>( thyra_P );

    // Create a composite operator C = (P + M + P^T)
    Teuchos::RCP<const Thyra::LinearOpBase<SC>> thyra_PpM =
        Thyra::add<SC>( thyra_P, thyra_M );
    Teuchos::RCP<const Thyra::LinearOpBase<SC>> thyra_C =
        Thyra::add<SC>( thyra_PpM, thyra_P_T );

    // Create parameters for stratimikos to setup the inverse operator.
    auto d_stratimikos_list = Teuchos::parameterList( "Stratimikos" );
    d_stratimikos_list->set( "Linear Solver Type", "Belos" );
    d_stratimikos_list->set( "Preconditioner Type", "None" );
    auto &linear_solver_types_list =
        d_stratimikos_list->sublist( "Linear Solver Types" );
    auto &belos_list = linear_solver_types_list.sublist( "Belos" );
    belos_list.set( "Solver Type", "Pseudo Block GMRES" );
    auto &solver_types_list = belos_list.sublist( "Solver Types" );
    auto &gmres_list = solver_types_list.sublist( "Pseudo Block GMRES" );
    gmres_list.set( "Convergence Tolerance", 1e-10 );
    gmres_list.set( "Verbosity",
                    Belos::Errors + Belos::Warnings + Belos::TimingDetails +
                        Belos::FinalSummary + Belos::StatusTestDetails );
    gmres_list.set( "Output Frequency", 1 );

    // Create the inverse of the composite operator C.
    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( d_stratimikos_list );
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC>> factory =
        Thyra::createLinearSolveStrategy( builder );
    Teuchos::RCP<const Thyra::LinearOpBase<SC>> thyra_C_inv =
        Thyra::inverse<SC>( *factory, thyra_C );

    // Create the composite operator B = (Q + N);
    Teuchos::RCP<const Thyra::LinearOpBase<SC>> thyra_B =
        Thyra::add<SC>( thyra_Q, thyra_N );

    // Create the coupling matrix A = (B * C^-1 * S).
    _thyra_operator = Thyra::multiply<SC>( thyra_B, thyra_C_inv, thyra_S );
    DTK_ENSURE( Teuchos::nonnull( _thyra_operator ) );

    _source = Teuchos::rcp( new Vector( S->getDomainMap(), 1 ) );
    _destination = Teuchos::rcp( new Vector( N->getRangeMap(), 1 ) );
    _thyra_X = Thyra::createMultiVector<SC>( _source );
    _thyra_Y = Thyra::createMultiVector<SC>( _destination );
}

template <typename DeviceType, typename CompactlySupportedRadialBasisFunction,
          typename PolynomialBasis>
void SplineOperator<DeviceType, CompactlySupportedRadialBasisFunction,
                    PolynomialBasis>::
    apply( Kokkos::View<double const *, DeviceType> source_values,
           Kokkos::View<double *, DeviceType> target_values ) const
{
    // Precondition: check that the source and the target are properly sized
    DTK_REQUIRE( source_values.extent( 0 ) ==
                 S->getDomainMap()->getNodeNumElements() );
    DTK_REQUIRE( target_values.extent( 0 ) ==
                 N->getRangeMap()->getNodeNumElements() );

    auto domain_map = S->getDomainMap();
    auto range_map = N->getRangeMap();

    auto source = Teuchos::rcp( new Vector( domain_map, 1 ) );
    auto destination = Teuchos::rcp( new Vector( range_map, 1 ) );

    Kokkos::deep_copy(
        Kokkos::subview( _source->getLocalViewDevice(), Kokkos::ALL, 0 ),
        source_values );

    _thyra_operator->apply( Thyra::NOTRANS, *_thyra_X, _thyra_Y.ptr(), 1, 0 );

    Kokkos::deep_copy(
        target_values,
        Kokkos::subview( _destination->getLocalViewDevice(), Kokkos::ALL, 0 ) );
}

} // end namespace DataTransferKit

// Explicit instantiation macro
#define DTK_SPLINE_OPERATOR_INSTANT( NODE )                                    \
    template class SplineOperator<typename NODE::device_type>;                 \
    template class SplineOperator<typename NODE::device_type, Wendland<0>,     \
                                  MultivariatePolynomialBasis<Quadratic, 3>>;

#endif
