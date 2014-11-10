//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_SpmdLinearOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief Map operator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAPOPERATOR_IMPL_HPP
#define DTK_MAPOPERATOR_IMPL_HPP

#include "DTK_DBC.hpp"

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Export.hpp>

#include <Thyra_SpmdMultiVectorBase.hpp>
#include <Thyra_SpmdVectorSpaceBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Scalar>
SpmdLinearOperator<Scalar>::SpmdLinearOperator(
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,int,std::size_t> >& crs_matrix )
    : d_crs_matrix( crs_matrix )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
template<class Scalar>
SpmdLinearOperator<Scalar>::~SpmdLinearOperator()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the range space.
template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
SpmdLinearOperator<Scalar>::range() const
{
    DTK_REQUIRE( Teuchos::nonnull(d_crs_matrix) );
    return Thyra::createVectorSpace( d_crs_matrix->getRangeMap() );
}

//---------------------------------------------------------------------------//
// Get the domain space.
template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
SpmdLinearOperator<Scalar>::domain() const
{
    DTK_REQUIRE( Teuchos::nonnull(d_crs_matrix) );
    return Thyra::createVectorSpace( d_crs_matrix->getDomainMap() );
}

//---------------------------------------------------------------------------//
// Clone the operator. The subclass must implement this function.
template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > 
SpmdLinearOperator<Scalar>::clone() const
{
    return Teuchos::rcp( new SpmdLinearOperator<Scalar>(d_crs_matrix) );
}
 
//---------------------------------------------------------------------------//
// Check if the given operator is supported.
template<class Scalar>
bool SpmdLinearOperator<Scalar>::opSupportedImpl( Thyra::EOpTransp M_trans ) const
{
    return ( M_trans == Thyra::NOTRANS );
}

//---------------------------------------------------------------------------//
// Apply the map operator to data defined on the entities by computing g =
// Minv*(v+A*f).
template<class Scalar>
void SpmdLinearOperator<Scalar>::applyImpl( 
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar>& X,
    const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> >& Y,
    const Scalar alpha,
    const Scalar beta ) const
{
    DTK_INSIST( opSupportedImpl(M_trans) );
    DTK_REQUIRE( Teuchos::nonnull(d_crs_matrix) );
    DTK_REQUIRE( Teuchos::nonnull(Y) );
  
    // Cast X and Y to SpmdMultiVectorBase objects in order to interrogate
    // their row decomposition.
    Teuchos::RCP<Thyra::SpmdMultiVectorBase<Scalar> > X_spmd =
	Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<Scalar> >(
	    Teuchos::rcpFromRef(X) );
    DTK_REQUIRE( Teuchos::nonnull(X_spmd) );
    Teuchos::RCP<Thyra::SpmdMultiVectorBase<Scalar> > Y_spmd =
	Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<Scalar> >(
	    Teuchos::rcpFromPtr(Y) );
    DTK_REQUIRE( Teuchos::nonnull(Y_spmd) );

    // Create a map for X.
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<Scalar> > X_space =
	X_spmd->spmdSpace();
    int num_X_elements = X_space->localSubDim();
    Teuchos::Array<std::size_t> X_elements( num_X_elements );
    for ( int i = 0; i < num_X_elements; ++i )
    {
	X_elements[i] = X_space->localOffset() + i;
    }
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > X_map =
	Tpetra::createNonContigMap<int,std::size_t>( 
	    X_elements, X_space->getComm() );

    // Create a map for Y.
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<Scalar> > Y_space =
	Y_spmd->spmdSpace();
    int num_Y_elements = Y_space->localSubDim();
    Teuchos::Array<std::size_t> Y_elements( num_Y_elements );
    for ( int i = 0; i < num_Y_elements; ++i )
    {
	Y_elements[i] = Y_space->localOffset() + i;
    }
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > Y_map =
	Tpetra::createNonContigMap<int,std::size_t>( 
	    Y_elements, Y_space->getComm() );

    // Create a copy of X.
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > X_tpetra =
	Tpetra::createMultiVector<Scalar,int,std::size_t>(
	    X_map, X_space->dim() );
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > X_copy =
	Thyra::createMultiVector( X_tpetra );
    Thyra::assign( X_copy.ptr(), X );

    // Create a copy of Y.
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > Y_tpetra =
	Tpetra::createMultiVector<Scalar,int,std::size_t>(
	    Y_map, Y_space->dim() );
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > Y_copy =
	Thyra::createMultiVector( Y_tpetra );
    Thyra::assign( Y_copy.ptr(), Y );

    // Export X to the operator domain map decomposition.
    Tpetra::Export X_export( X_map, d_crs_matrix->getDomainMap() );
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > X_domain_dist =
	Tpetra::createMultiVector<Scalar,int,std::size_t>(
	    d_crsMatrix->getDomainMap(), X_space->dim() );
    X_domain_dist->doExport( *X_tpetra, X_export, Tpetra::INSERT );

    // Export Y to the operator range decomposition.
    Tpetra::Export Y_export( Y_map, d_crs_matrix->getRangeMap() );
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > Y_range_dist =
	Tpetra::createMultiVector<Scalar,int,std::size_t>(
	    d_crsMatrix->getRangeMap(), Y_space->dim() );
    Y_range_dist->doExport( *Y_tpetra, Y_export, Tpetra::INSERT );

    // Apply the operator.
    d_crs_matrix->apply( *X_domain_dist, *Y_range_dist, 
			 Teuchos::NO_TRANS, alpha, beta );

    // Import Y to the original decomposition.
    Y_tpetra.doImport( *Y_range_dist, Y_export, Tpetra::INSERT );
    Thyra::assign( Y, *Y_copy );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

# endif // end DTK_MAPOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SpmdLinearOperator_impl.hpp
//---------------------------------------------------------------------------//
