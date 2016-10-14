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
 * \brief DTK_MapOperator.hpp
 * \author Stuart R. Slattery
 * \brief Map operator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAPOPERATOR_HPP
#define DTK_MAPOPERATOR_HPP

#include "DTK_FunctionSpace.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MapOperator
  \brief Map operator interface.

  A map operator maps a field in one entity set to another entity set.
*/
//---------------------------------------------------------------------------//
class MapOperator : public Tpetra::Operator<double, int, SupportId>
{
  public:
    //! Root class typedef.
    typedef Tpetra::Operator<double, int, SupportId> Root;

    //! Map typedef.
    typedef Tpetra::Map<int, SupportId, typename Root::node_type> TpetraMap;

    //! MultiVector typedef.
    typedef Tpetra::MultiVector<double, int, SupportId,
                                typename Root::node_type>
        TpetraMultiVector;

    /*!
     * \brief Constructor.
     *
     * \param domain_map Parallel map for domain vectors this map should be
     * compatible with.
     *
     * \param range_map Parallel map for range vectors this map should be
     * compatible with.
     */
    MapOperator( const Teuchos::RCP<const TpetraMap> &domain_map,
                 const Teuchos::RCP<const TpetraMap> &range_map );

    /*!
     * \brief Destructor.
     */
    virtual ~MapOperator();

    /*
     * \brief Setup the map operator from a domain function space and a range
     * function space. Subclasses should override the setupImpl() function.
     *
     * \param domain_function The function that contains the data that will be
     * sent to the range. Must always be nonnull but the pointers it contains
     * may be null of no entities are on-process.
     *
     * \param range_space The function that will receive the data from the
     * domain. Must always be nonnull but the pointers it contains
     * may be null of no entities are on-process.
     *
     * \param parameters Parameters for the setup.
     */
    void setup( const Teuchos::RCP<FunctionSpace> &domain_space,
                const Teuchos::RCP<FunctionSpace> &range_space );

    /*!
     * \brief Return whether or not the operator has been setup.
     */
    bool setupIsComplete() const;

    //@{
    //! Tpetra::Operator interface.
    Teuchos::RCP<const TpetraMap> getDomainMap() const override;
    Teuchos::RCP<const TpetraMap> getRangeMap() const override;
    void
    apply( const TpetraMultiVector &X, TpetraMultiVector &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           double alpha = Teuchos::ScalarTraits<double>::one(),
           double beta = Teuchos::ScalarTraits<double>::zero() ) const override;
    bool hasTransposeApply() const override;
    //@}

  protected:
    //! Tranpose apply option.
    virtual bool hasTransposeApplyImpl() const = 0;

    //! Setup implementation. Subclasses should override.
    virtual void
    setupImpl( const Teuchos::RCP<FunctionSpace> &domain_space,
               const Teuchos::RCP<FunctionSpace> &range_space ) = 0;

    //! Apply implementation. Subclasses should override.
    virtual void applyImpl( const TpetraMultiVector &X, TpetraMultiVector &Y,
                            Teuchos::ETransp mode, double alpha,
                            double beta ) const = 0;

  private:
    //! Domain map.
    Teuchos::RCP<const TpetraMap> d_domain_map;

    //! Range map.
    Teuchos::RCP<const TpetraMap> d_range_map;

    //! True if setup has been completed.
    bool d_setup_is_complete;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MAPOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MapOperator.hpp
//---------------------------------------------------------------------------//
