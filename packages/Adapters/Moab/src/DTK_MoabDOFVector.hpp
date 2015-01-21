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
 * \brief DTK_MoabDOFVector.hpp
 * \author Stuart R. Slattery
 * \brief Helper functions for managing Moab mesh DOF vectors.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABDOFVECTOR_HPP
#define DTK_MOABDOFVECTOR_HPP

#include <MBParallelComm.hpp>

#include <Teuchos_RCP.hpp>

#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabDOFVector
  \brief Helper functions for managing Moab mesh DOF vectors.
*/
//---------------------------------------------------------------------------//
class MoabDOFVector
{
  public:

    /*!
     * \brief Constructor.
     */
    MoabDOFVector()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~MoabDOFVector()
    { /* ... */ }

    /*!
     * \brief Given a Moab tag, create a Tpetra vector that maps to the tag
     * DOFs on the given mesh set and pull the data from the tag.
     * \param moab_mesh The mesh over which the tag is defined.
     * \param mesh_set The mesh set over which the tag is defined.
     * \param tag The moab tag.
     * \return A Tpetra MultiVector indexed according to the tag entities
     * with a vector for each tag dimension.
     */
    template<class Scalar>
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
    pullTpetraMultiVectorFromMoabTag( const moab::ParallelComm& moab_mesh,
				      const moab::EntityHandle& mesh_set,
				      const moab::Tag& tag );

    /*!
     * \brief Given a Tpetra vector of DOF data, push the data into a given
     * Moab tag on the given mesh set.
     * \param tag_dofs A Tpetra vector containing the tag DOFs. One vector
     * for each dimension of the tag.
     * \param moab_mesh The mesh over which the tag is defined.
     * \param mesh_set The mesh set over which the tag is defined.
     * \param tag The Moab tag.
     */
    template<class Scalar>
    static void pushTpetraMultiVectorToMoabTag(
	const Tpetra::MultiVector<Scalar,int,std::size_t>& tag_dofs,
	const moab::ParallelComm& moab_mesh,
	const moab::EntityHandle& mesh_set,
	const moab::Tag& tag );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_MoabDOFVector_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_MOABDOFVECTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabDOFVector.hpp
//---------------------------------------------------------------------------//
