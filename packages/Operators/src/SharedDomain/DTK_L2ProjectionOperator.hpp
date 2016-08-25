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
 * \brief DTK_L2ProjectionOperator.hpp
 * \author Stuart R. Slattery
 * \brief L2 projection operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_L2PROJECTIONOPERATOR_HPP
#define DTK_L2PROJECTIONOPERATOR_HPP

#include "DTK_MapOperator.hpp"
#include "DTK_Types.hpp"
#include "DTK_EntityIterator.hpp"
#include "DTK_IntegrationPointSet.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_CrsMatrix.hpp>

#include <Thyra_LinearOpBase.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class L2ProjectionOperator
  \brief L2 projection operator.

  Constructs and solves the L2 projection problem on a shared domain. The
  Galerkin problem is assembled over the range (target) entity set.
*/
//---------------------------------------------------------------------------//
class L2ProjectionOperator : virtual public MapOperator
{
  public:

    //! Root class tyepdef.
    typedef MapOperator Base;
    typedef typename Base::Root Root;
    typedef typename Root::scalar_type Scalar;
    typedef typename Root::local_ordinal_type LO;
    typedef typename Root::global_ordinal_type GO;
    typedef typename Base::TpetraMultiVector TpetraMultiVector;
    typedef typename Base::TpetraMap TpetraMap;
    
    /*!
     * \brief Constructor.
     */
    L2ProjectionOperator(
	const Teuchos::RCP<const TpetraMap>& domain_map,
	const Teuchos::RCP<const TpetraMap>& range_map,
	const Teuchos::ParameterList& parameters );

  protected:

    /*
     * \brief Setup the map operator from a domain entity set and a range
     * entity set.
     *
     * \param domain_map Parallel map for domain vectors this map should be
     * compatible with.
     *
     * \param domain_function The function that contains the data that will be
     * sent to the range. Must always be nonnull but the pointers it contains
     * may be null of no entities are on-process.
     *
     * \param range_map Parallel map for range vectors this map should be
     * compatible with.
     *
     * \param range_space The function that will receive the data from the
     * domain. Must always be nonnull but the pointers it contains to entity
     * data may be null of no entities are on-process.
     *
     * \param parameters Parameters for the setup.
     */
    void setupImpl( const Teuchos::RCP<FunctionSpace>& domain_space,
		    const Teuchos::RCP<FunctionSpace>& range_space ) override;

    /*!
     * \brief Apply the operator.
     */
    void applyImpl(
	const TpetraMultiVector& X,
	TpetraMultiVector &Y,
	Teuchos::ETransp mode = Teuchos::NO_TRANS,
	double alpha = Teuchos::ScalarTraits<double>::one(),
	double beta = Teuchos::ScalarTraits<double>::zero()) const override;

    /*
     * \brief Transpose apply option.
     */
    bool hasTransposeApplyImpl() const override;

  private:

    // Assemble the mass matrix and range integration point set.
    void assembleMassMatrix(
	const Teuchos::RCP<FunctionSpace>& range_space,
	EntityIterator range_iterator,
	Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO> >& mass_matrix,
	Teuchos::RCP<IntegrationPointSet>& range_ip_set );

    // Assemble the coupling matrix.
    void assembleCouplingMatrix(
	const Teuchos::RCP<FunctionSpace>& domain_space,
	EntityIterator domain_iterator,
	const Teuchos::RCP<IntegrationPointSet>& range_ip_set,
	Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO> >& coupling_matrix );
    
  private:

    // Order of numerical integration for assembly of the Galerkin problem.
    int d_int_order;
    
    // Search sublist.
    Teuchos::ParameterList d_search_list;

    // Coupling matrix.
    Teuchos::RCP<const Thyra::LinearOpBase<double> > d_l2_operator;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_L2PROJECTIONOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_L2ProjectionOperator.hpp
//---------------------------------------------------------------------------//
