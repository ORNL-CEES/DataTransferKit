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
 * \file   DTK_MovingLeastSquareReconstructionOperator.hpp
 * \author Stuart R. Slattery
 * \brief Parallel moving least square interpolator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOVINGLEASTSQUARERECONSTRUCTIONOPERATOR_HPP
#define DTK_MOVINGLEASTSQUARERECONSTRUCTIONOPERATOR_HPP

#include "DTK_MapOperator.hpp"
#include "DTK_RadialBasisPolicy.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_CrsMatrix.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class MovingLeastSquareReconstructionOperator
 * \brief Parallel moving least square interpolator MapOperator
 * implementation.
 */
//---------------------------------------------------------------------------//
template<class Basis,int DIM>
class MovingLeastSquareReconstructionOperator : virtual public MapOperator
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

    /*!
     * \brief Constructor.
     *
     * \param domain_map Parallel map for domain vectors this map should be
     * compatible with.
     *
     * \param range_map Parallel map for range vectors this map should be
     * compatible with.
     */
    MovingLeastSquareReconstructionOperator(
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

    // Extract node coordinates and ids from an iterator.
    void getNodeCoordsAndIds( const Teuchos::RCP<FunctionSpace>& space,
                              EntityIterator iterator,
                              Teuchos::ArrayRCP<double>& centers,
                              Teuchos::ArrayRCP<GO>& support_ids ) const;
    
  private:

    // Flag for search type. True if kNN, false if radius.
    bool d_use_knn;

    // k-nearest-neighbors for support.
    int d_knn;
    
    // Basis radius for support.
    double d_radius;

    // Domain entity topological dimension. Default is 0 (vertex).
    int d_domain_entity_dim;

    // Range entity topological dimension. Default is 0 (vertex).
    int d_range_entity_dim;

    // Coupling matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> > d_coupling_matrix;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOVINGLEASTSQUARERECONSTRUCTIONOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MovingLeastSquareReconstructionOperator.hpp
//---------------------------------------------------------------------------//

