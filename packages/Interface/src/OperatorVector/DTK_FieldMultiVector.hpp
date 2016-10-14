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
 * \brief DTK_FieldMultiVector.hpp
 * \author Stuart R. Slattery
 * \brief MultiVector interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDMULTIVECTOR_HPP
#define DTK_FIELDMULTIVECTOR_HPP

#include "DTK_EntitySet.hpp"
#include "DTK_Field.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class FieldMultiVector
  \brief MultiVector interface.

  FieldMultiVector provides a Tpetra::MultiVector wrapper around application
  field data. Client implementations of the Field interface provide read/write
  access to field data on an entity-by-entity basis. The FieldMultiVector then
  manages the copying of data between the application and the Tpetra vector
  using the client implementations for data access.
*/
//---------------------------------------------------------------------------//
class FieldMultiVector : public Tpetra::MultiVector<double, int, SupportId>
{
  public:
    //! MultiVector typedef.
    typedef Tpetra::MultiVector<double, int, SupportId> Base;
    typedef typename Base::local_ordinal_type LO;
    typedef typename Base::global_ordinal_type GO;

    /*!
     * \brief Comm constructor. This will allocate the Tpetra vector.
     *
     * \param field The field for which we are building a vector.
     *
     * \param global_comm The global communicator over which the field is
     * defined.
     */
    FieldMultiVector( const Teuchos::RCP<const Teuchos::Comm<int>> &global_comm,
                      const Teuchos::RCP<Field> &field );

    /*!
     * \brief Entity set constructor. This will allocate the Tpetra vector.
     *
     * \param field The field for which we are building a vector.
     *
     * \param entity_set The entity set over which the field is defined.
     */
    FieldMultiVector( const Teuchos::RCP<Field> &field,
                      const Teuchos::RCP<const EntitySet> &entity_set );

    /*!
     * \brief Pull data from the application and put it in the vector.
     */
    void pullDataFromApplication();

    /*!
     * \brief Push data from the vector into the application.
     */
    void pushDataToApplication();

  private:
    // The field this multivector is managing.
    Teuchos::RCP<Field> d_field;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_FIELDMULTIVECTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldMultiVector.hpp
//---------------------------------------------------------------------------//
