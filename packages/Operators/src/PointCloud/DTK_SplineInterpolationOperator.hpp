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
 * \file   DTK_SplineInterpolationOperator.hpp
 * \author Stuart R. Slattery
 * \brief Parallel spline interpolator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEINTERPOLATIONOPERATOR_HPP
#define DTK_SPLINEINTERPOLATIONOPERATOR_HPP

#include "DTK_MapOperator.hpp"
#include "DTK_RadialBasisPolicy.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_Map.hpp>

#include <Thyra_LinearOpBase.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class SplineInterpolationOperator
 * \brief Parallel spline interpolator.
 *
 * The SplineInterpolationOperator is the top-level driver for parallel interpolation
 * problems.
 */
//---------------------------------------------------------------------------//
template<class Basis,int DIM>
class SplineInterpolationOperator : virtual public MapOperator
{
  public:

    //@{
    //! Typedefs.
    typedef MapOperator Base;
    typedef typename Base::Root Root;
    typedef typename Root::scalar_type Scalar;
    typedef typename Root::local_ordinal_type LO;
    typedef typename Root::global_ordinal_type GO;
    typedef typename Base::TpetraMultiVector TpetraMultiVector;
    typedef typename Base::TpetraMap TpetraMap;
    typedef RadialBasisPolicy<Basis> BP;
    //@}

    /*
     * \brief Constructor.
     *
     * \param domain_map Parallel map for domain vectors this map should be
     * compatible with.
     *
     * \param range_map Parallel map for range vectors this map should be
     * compatible with.
     */
    SplineInterpolationOperator(     
	const Teuchos::RCP<const TpetraMap>& domain_map,
	const Teuchos::RCP<const TpetraMap>& range_map,
	const Teuchos::ParameterList& parameters );

  protected:
    
    /*
     * \brief Setup the map operator from a domain entity set and a range
     * entity set.
     *
     * \param domain_function The function that contains the data that will be
     * sent to the range. Must always be nonnull but the pointers it contains
     * may be null of no entities are on-process.
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
    
  private:

    // Build the concrete operators.
    void buildConcreteOperators(
	const Teuchos::RCP<FunctionSpace>& domain_space,
	const Teuchos::RCP<FunctionSpace>& range_space,
	Teuchos::RCP<const Root>& S,
	Teuchos::RCP<const Root>& P,
	Teuchos::RCP<const Root>& M,
	Teuchos::RCP<const Root>& Q,
	Teuchos::RCP<const Root>& N ) const;

  private:

    // Basis radius.
    double d_radius;

    // Coupling matrix.
    Teuchos::RCP<const Thyra::LinearOpBase<double> > d_coupling_matrix;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_SplineInterpolationOperator_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEINTERPOLATIONOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineInterpolationOperator.hpp
//---------------------------------------------------------------------------//

