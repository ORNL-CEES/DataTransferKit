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
 * \brief DTK_EntityCenteredDOFVector_impl.hpp
 * \author Stuart R. Slattery
 * \brief Entity-centered DOF vector.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYCENTEREDDOFVECTOR_IMPL_HPP
#define DTK_ENTITYCENTEREDDOFVECTOR_IMPL_HPP

#include "DTK_DBC.hpp"

#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Map.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template<class Scalar>
Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::ArrayView<Entity>& entities,
    const int field_dim,
    const Teuchos::ArrayView<Scalar>& dof_data )
{
    std::size_t lda = entities.size();
    DTK_REQUIRE( lda*field_dim == Teuchos::as<std::size_t>(dof_data.size()) );

    // Extract the entity ids.
    Teuchos::Array<std::size_t> entity_ids( lda );
    auto id_it = entity_ids.begin();
    for ( auto entity_it = entities.begin();
	  entity_it != entities.end();
	  ++entity_it, ++id_it )
    {
	*id_it = entity_it->id();
    }
	  
    // Construct a map.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	Tpetra::createNonContigMap<int,std::size_t>( entity_ids(), comm );

    // Build a tpetra multivector.
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > vector =
	Tpetra::createMultiVector<Scalar,int,std::size_t>( map, field_dim );

    // Copy the data.
    Teuchos::ArrayRCP<Scalar> vector_view = vector->get1dViewNonConst();
    std::copy( dof_data.begin(), dof_data.end(), vector_view.begin() );
    return vector;
}

//---------------------------------------------------------------------------//
template<class Scalar>
void EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
    const Tpetra::MultiVector<Scalar,int,std::size_t>& vector,
    Teuchos::ArrayView<Scalar>&& dof_data )
{
    Teuchos::ArrayRCP<const Scalar> vector_view = vector.get1dView();
    std::copy( vector_view.begin(), vector_view.end(), dof_data.begin() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYCENTEREDDOFVECTOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityCenteredDOFVector_impl.hpp
//---------------------------------------------------------------------------//
