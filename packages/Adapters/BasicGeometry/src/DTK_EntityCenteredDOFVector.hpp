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
 * \brief DTK_EntityCenteredDOFVector.hpp
 * \author Stuart R. Slattery
 * \brief Entity-centered DOF vector.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYCENTEREDDOFVECTOR_HPP
#define DTK_ENTITYCENTEREDDOFVECTOR_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>

#include <Thyra_MultiVectorBase.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityCenteredDOFVector
  \brief Entity-centered DOF vector.

  Tools for constructing vectors for data bound to entity-centered shape
  functions. 
*/
//---------------------------------------------------------------------------//
class EntityCenteredDOFVector
{
  public:

    /*!
     * \brief Constructor.
     */
    EntityCenteredDOFVector()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~EntityCenteredDOFVector()
    { /* ... */ }

    /*!
     * \brief Given a set of entity ids and DOF data bound to the center of
     * those entites, build a Thyra vector.
     * \param comm The communicator the DOFs are defined over.
     * \param entity_ids The ids of the entities the DOFs are defined over
     * (cast as std::size_t instead of EntityId).
     * \param dof_data Reference counted array of the DOF data for the given
     * entity ids.
     * \param lda Single vecto length. Should be the same as the size of
     * entity_ids.
     * \param num_vectors The number of vectors in the multivector.
     */
    template<class Scalar>
    static Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > 
    createThyraMultiVector( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
			    const Teuchos::ArrayView<const std::size_t>& entity_ids,
			    const Teuchos::ArrayRCP<Scalar>& dof_data,
			    const std::size_t lda,
			    const std::size_t num_vectors );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYCENTEREDDOFVECTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityCenteredDOFVector.hpp
//---------------------------------------------------------------------------//
