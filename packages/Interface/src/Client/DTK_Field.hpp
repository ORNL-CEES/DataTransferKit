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
 * \brief DTK_Field.hpp
 * \author Stuart R. Slattery
 * \brief Application field interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELD_HPP
#define DTK_FIELD_HPP

#include "DTK_Types.hpp"

#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class Field
  \brief Field interface.

  Field provides an access wrapper around application field data. Client API
  implementations provided read/write access to field data on an
  entity-by-entity basis. Field is a translation layer between entity support
  ids provided by shape functions and the actual indexing into application
  data for a given quantity of interest. For example, multiple nodal fields
  will have entity support ids that are identical while the data in the
  underlying application will be indexed differently.

  We need this extra layer of indirection because MapOperator may be
  constructed to be compatible with multiple subsets of a single vector in an
  application. For example, if nodal data of pressure and temperature were in
  a single state vector in an application, two Field objects could be
  constructed, one for each field. A single nodal MapOperator could then be
  applied to each vector constructed from the fields to independently perform
  the solution transfer.
*/
//---------------------------------------------------------------------------//
class Field
{
  public:
    /*!
     * \brief Constructor.
     */
    Field() { /* ... */}

    /*!
     * \brief Destructor.
     */
    virtual ~Field() { /* ... */}

    /*!
     * \brief Get the dimension of the field.
     */
    virtual int dimension() const = 0;

    /*!
     * \brief Get the locally-owned support location ids of the field.
     */
    virtual Teuchos::ArrayView<const SupportId> getLocalSupportIds() const = 0;

    /*!
     * \brief Given a local support id and a dimension, read data from the
     * application field.
     */
    virtual double readFieldData( const SupportId support_id,
                                  const int dimension ) const = 0;

    /*!
     * \brief Given a local support id, dimension, and field value, write data
     * into the application field.
     */
    virtual void writeFieldData( const SupportId support_id,
                                 const int dimension, const double data ) = 0;

    /*!
     * \brief Finalize a field after writing into it. This lets some clients
     * do a post-process (e.g. update ghost values). Default finalize does
     * nothing.
     */
    virtual void finalizeAfterWrite() { /* ... */}
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_FIELD_HPP

//---------------------------------------------------------------------------//
// end DTK_Field.hpp
//---------------------------------------------------------------------------//
