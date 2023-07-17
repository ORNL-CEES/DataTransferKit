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

#ifndef DTK_SPLINE_PROLONGATION_OPERATOR_HPP
#define DTK_SPLINE_PROLONGATION_OPERATOR_HPP

#include <DTK_Types.h>

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

#include <Trilinos_version.h>

namespace DataTransferKit
{

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal,
          typename Node>
class SplineProlongationOperator
    : public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
    using Map = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using MultiVector =
        Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  public:
    // Constructor.
    SplineProlongationOperator( const int num_polynomial_dofs,
                                const Teuchos::RCP<const Map> &domain_map )
        : _domain_map( domain_map )
    {
        // Create a range map.
        Teuchos::ArrayView<const GlobalOrdinal> domain_elements =
            _domain_map->getLocalElementList();
        _lda = domain_elements.size();

        const auto old_size = domain_elements.size();
        Teuchos::Array<GlobalOrdinal> global_ids( old_size +
                                                  num_polynomial_dofs );
        if ( _domain_map->getComm()->getRank() == 0 )
        {
            GlobalOrdinal max_id = _domain_map->getMaxAllGlobalIndex() + 1;

            global_ids( 0, old_size ).assign( domain_elements );
            for ( int i = 0; i < num_polynomial_dofs; ++i )
            {
                global_ids[old_size + i] = max_id + i;
            }
            domain_elements = global_ids();
        }
        _range_map = Tpetra::createNonContigMapWithNode<LocalOrdinal,
                                                        GlobalOrdinal, Node>(
            domain_elements, _domain_map->getComm() );
        DTK_ENSURE( Teuchos::nonnull( _range_map ) );
    }

    //! The Map associated with the domain of this operator, which must be
    //! compatible with X.getMap().
    Teuchos::RCP<const Map> getDomainMap() const override
    {
        return _domain_map;
    }

    //! The Map associated with the range of this operator, which must be
    //! compatible with Y.getMap().
    Teuchos::RCP<const Map> getRangeMap() const override { return _range_map; }

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X +
        \beta \cdot Y\f$. However, the details of operation vary according to
        the values of \c alpha and \c beta. Specifically - if <tt>beta ==
        0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y
        (including NaNs) are ignored.  - if <tt>alpha == 0</tt>, apply()
        <b>may</b> short-circuit the operator, so that any values in \c X
        (including NaNs) are ignored.
     */
    void
    apply( const MultiVector &X, MultiVector &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero() ) const override
    {
        DTK_REQUIRE( _domain_map->isSameAs( *( X.getMap() ) ) );
        DTK_REQUIRE( _range_map->isSameAs( *( Y.getMap() ) ) );
        DTK_REQUIRE( X.getNumVectors() == Y.getNumVectors() );

        using DeviceType = typename Node::device_type;
        using ExecutionSpace = typename DeviceType::execution_space;

        auto x_view = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
        auto y_view = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);

        auto const num_vectors = x_view.extent_int( 1 );

        Y.scale( beta );
        Kokkos::parallel_for( DTK_MARK_REGION( "spline_prolongation::apply" ),
                              Kokkos::RangePolicy<ExecutionSpace>( 0, _lda ),
                              KOKKOS_LAMBDA( int const i ) {
                                  for ( int j = 0; j < num_vectors; ++j )
                                      y_view( i, j ) += alpha * x_view( i, j );
                              } );
    }
    /// \brief Whether this operator supports applying the transpose or
    /// conjugate transpose.
    bool hasTransposeApply() const override { return false; }

  private:
    int _lda;
    Teuchos::RCP<const Map> _domain_map;
    Teuchos::RCP<const Map> _range_map;
};

} // end namespace DataTransferKit

#endif
