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

#ifndef DTK_POLYNOMIAL_MATRIX_HPP
#define DTK_POLYNOMIAL_MATRIX_HPP

#include <DTK_Types.h>

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <DTK_DBC.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <mpi.h>
#endif

#include <Tpetra_Export.hpp>

#include <Trilinos_version.h>

namespace DataTransferKit
{

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal,
          typename Node>
class PolynomialMatrix
    : public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
    using DeviceType = typename Node::device_type;
    using Map = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using MultiVector =
        Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  public:
    PolynomialMatrix( Kokkos::View<double **, DeviceType> vandermonde,
                      const Teuchos::RCP<const Map> &domain_map,
                      const Teuchos::RCP<const Map> &range_map )
        : _vandermonde( vandermonde )
        , _domain_map( domain_map )
        , _range_map( range_map )
    {
    }

    Teuchos::RCP<const Map> getDomainMap() const override
    {
        return _domain_map;
    }

    Teuchos::RCP<const Map> getRangeMap() const override { return _range_map; }

    void
    apply( const MultiVector &X, MultiVector &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero() ) const override
    {
        DTK_REQUIRE( _domain_map->isSameAs( *( X.getMap() ) ) );
        DTK_REQUIRE( _range_map->isSameAs( *( Y.getMap() ) ) );
        DTK_REQUIRE( X.getNumVectors() == Y.getNumVectors() );

        using ExecutionSpace = typename DeviceType::execution_space;

        auto comm = _domain_map->getComm();
#ifdef HAVE_MPI
        Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
            Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( comm );
        Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
            mpi_comm->getRawMpiComm();
        MPI_Comm raw_comm = ( *opaque_comm )();
#endif

        // Get the size of the problem and view of the local vectors.
        int const local_length = _vandermonde.extent( 0 );
        int const poly_size = _vandermonde.extent( 1 );
        int const num_vec = X.getNumVectors();

        // To avoid capturing *this
        auto vandermonde = _vandermonde;

        Y.scale( beta );

        if ( mode == Teuchos::NO_TRANS )
        {
            Kokkos::View<double **, DeviceType> x_poly( "x_poly", poly_size,
                                                        num_vec );
            if ( 0 == comm()->getRank() )
            {
                auto x_view = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
                auto const n = x_view.extent( 0 );
                Kokkos::deep_copy(
                    x_poly, Kokkos::subview(
                                x_view, Kokkos::make_pair( n - poly_size, n ),
                                Kokkos::ALL ) );
            }

#ifdef HAVE_MPI
            {
                // Broadcast the polynomial components of X from the root rank.
                auto x_poly_host = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace{}, x_poly );
                MPI_Bcast( x_poly_host.data(), poly_size * num_vec, MPI_DOUBLE,
                           0, raw_comm );
                Kokkos::deep_copy( x_poly, x_poly_host );
            }
#endif
            auto y_view = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);
            Kokkos::parallel_for(
                DTK_MARK_REGION( "polynomial_matrix::apply::no_trans" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, local_length ),
                KOKKOS_LAMBDA( int const i ) {
                    for ( int j = 0; j < num_vec; ++j )
                        for ( int p = 0; p < poly_size; ++p )
                            y_view( i, j ) +=
                                alpha * vandermonde( i, p ) * x_poly( p, j );
                } );
        }
        else if ( mode == Teuchos::TRANS )
        {
            // Export X to the polynomial decomposition.
            MultiVector work( _range_map, num_vec );
            Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node> exporter(
                _domain_map, _range_map );
            work.doExport( X, exporter, Tpetra::INSERT );

            // Do the local mat-vec.
            auto work_view = work.getLocalViewDevice(Tpetra::Access::ReadOnly);
            Kokkos::View<double **, DeviceType> products( "products", poly_size,
                                                          num_vec );
            {
                auto scatter_products =
                    Kokkos::Experimental::create_scatter_view( products );
                Kokkos::parallel_for(
                    DTK_MARK_REGION( "polynomial_matrix::apply::trans" ),
                    Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<3>>(
                        {0, 0, 0}, {local_length, num_vec, poly_size} ),
                    KOKKOS_LAMBDA( int const i, int const j, int const p ) {
                        auto access = scatter_products.access();
                        access( p, j ) +=
                            alpha * vandermonde( i, p ) * work_view( i, j );
                    } );
                Kokkos::Experimental::contribute( products, scatter_products );
            }

            // Reduce the results back to the root rank.
            Kokkos::View<double **, DeviceType> product_sums(
                "product_sums", poly_size, num_vec );
#ifdef HAVE_MPI
            {
                auto products_host = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace{}, products );
                auto product_sums_host = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace{}, product_sums );
                MPI_Reduce( products_host.data(), product_sums_host.data(),
                            poly_size * num_vec, MPI_DOUBLE, MPI_SUM, 0,
                            raw_comm );
                Kokkos::deep_copy( product_sums, product_sums_host );
            }
#else
            product_sums = products;
#endif

            // Assign the values to Y on the root rank.
            // Note: no alpha here as we used it above.
            if ( 0 == comm->getRank() )
            {
                auto y_view = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);

                auto const n = y_view.extent( 0 );
                Kokkos::deep_copy(
                    Kokkos::subview( y_view,
                                     Kokkos::make_pair( n - poly_size, n ),
                                     Kokkos::ALL ),
                    product_sums );
            }
        }
    }

    bool hasTransposeApply() const override { return true; }

  private:
    Kokkos::View<double **, DeviceType> _vandermonde;

    Teuchos::RCP<const Map> _domain_map;
    Teuchos::RCP<const Map> _range_map;
};

} // end namespace DataTransferKit

#endif
