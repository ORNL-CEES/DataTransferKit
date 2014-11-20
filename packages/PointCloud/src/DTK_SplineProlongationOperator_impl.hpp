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
 * \file   DTK_SplineProlongationOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Spline transformation operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEPROLONGATIONOPERATOR_IMPL_HPP
#define DTK_SPLINEPROLONGATIONOPERATOR_IMPL_HPP

#include "DTK_DBC.hpp"

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar,class GO>
SplineProlongationOperator<Scalar,GO>::SplineProlongationOperator(
    const int offset,
    const Teuchos::RCP<const Tpetra::Map<int,GO> >& domain_map )
    : d_offset( offset )
    , d_domain_map( domain_map )
{
    // Create a range map.
    Teuchos::ArrayView<const GO> domain_elements = 
	d_domain_map->getNodeElementList();
    d_lda = domain_elements.size();
    Teuchos::Array<GO> global_ids;
    if ( d_domain_map->getComm()->getRank() == 0 )
    {
	global_ids.resize( d_offset + domain_elements.size() );
	global_ids( d_offset, domain_elements.size() ).assign( domain_elements );
	GO max_id = d_domain_map->getMaxAllGlobalIndex() + 1;
	for ( int i = 0; i < d_offset; ++i )
	{
	    global_ids[i] = max_id + i;
	}
	domain_elements = global_ids();
    }
    else
    {
	d_offset = 0;
    }
    d_range_map = Tpetra::createNonContigMap<int,GO>(
	domain_elements, d_domain_map->getComm() );
    DTK_ENSURE( Teuchos::nonnull(d_range_map) );
}

//---------------------------------------------------------------------------//
// Apply operation. 
template<class Scalar,class GO>
void SplineProlongationOperator<Scalar,GO>::apply(
    const Tpetra::MultiVector<Scalar,int,GO> &X,
    Tpetra::MultiVector<Scalar,int,GO> &Y,
    Teuchos::ETransp mode,
    Scalar alpha,
    Scalar beta ) const
{
    DTK_REQUIRE( d_domain_map->isSameAs(*(X.getMap())) );
    DTK_REQUIRE( d_range_map->isSameAs(*(Y.getMap())) );
    DTK_REQUIRE( X.getNumVectors() == Y.getNumVectors() );

    Y.putScalar( 0.0 );
    
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > X_view = X.get2dView();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > Y_view = Y.get2dViewNonConst();
    for ( int n = 0; n < X.getNumVectors(); ++n )
    {
	Y_view[n](d_offset,d_lda).assign( X_view[n]() );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEPROLONGATIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineProlongationOperator_impl.hpp
//---------------------------------------------------------------------------//

