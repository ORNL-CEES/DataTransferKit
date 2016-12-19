//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file   DTK_SplineProlongationOperator.hpp
 * \author Stuart R. Slattery
 * \brief  Spline transformation operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEPROLONGATIONOPERATOR_HPP
#define DTK_SPLINEPROLONGATIONOPERATOR_HPP

#include "DTK_MapOperator.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class SplineProlongationOperator
 * \brief Prolongation operator for projecting a vector into the extended
 * spline space.
 */
//---------------------------------------------------------------------------//
class SplineProlongationOperator : public MapOperator::Root
{
  public:
    //@{
    //! Typedefs.
    typedef typename MapOperator::Scalar Scalar;
    typedef typename MapOperator::LO LO;
    typedef typename MapOperator::GO GO;
    typedef typename MapOperator::Node Node;
    //@}

    // Constructor.
    SplineProlongationOperator(
        const int offset,
        const Teuchos::RCP<const Tpetra::Map<LO, GO, Node>> &domain_map );

    //! The Map associated with the domain of this operator, which must be
    //! compatible with X.getMap().
    Teuchos::RCP<const Tpetra::Map<LO, GO, Node>> getDomainMap() const override
    {
        return d_domain_map;
    }

    //! The Map associated with the range of this operator, which must be
    //! compatible with Y.getMap().
    Teuchos::RCP<const Tpetra::Map<LO, GO, Node>> getRangeMap() const override
    {
        return d_range_map;
    }

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
    apply( const Tpetra::MultiVector<Scalar, LO, GO, Node> &X,
           Tpetra::MultiVector<Scalar, LO, GO, Node> &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero() ) const override;

    /// \brief Whether this operator supports applying the transpose or
    /// conjugate transpose.
    bool hasTransposeApply() const override { return false; }

  private:
    // Prolongation offset.
    int d_offset;

    // LDA
    int d_lda;

    // Domain map.
    Teuchos::RCP<const Tpetra::Map<LO, GO, Node>> d_domain_map;

    // Range map.
    Teuchos::RCP<const Tpetra::Map<LO, GO, Node>> d_range_map;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEPROLONGATIONOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineProlongationOperator.hpp
//---------------------------------------------------------------------------//
