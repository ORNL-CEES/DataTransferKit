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
 * \brief DTK_FieldMultiVector.cpp
 * \author Stuart R. Slattery
 * \brief MultiVector interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_FieldMultiVector.hpp"
#include "DTK_DBC.hpp"

#include <Tpetra_Map.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Comm constructor.
FieldMultiVector::FieldMultiVector(
    const Teuchos::RCP<const Teuchos::Comm<int>> &global_comm,
    const Teuchos::RCP<Field> &field )
    : Base( Tpetra::createNonContigMap<int, SupportId>(
                field->getLocalSupportIds(), global_comm ),
            field->dimension() )
    , d_field( field )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Entity set constructor.
FieldMultiVector::FieldMultiVector(
    const Teuchos::RCP<Field> &field,
    const Teuchos::RCP<const EntitySet> &entity_set )
    : Base( Tpetra::createNonContigMap<int, SupportId>(
                field->getLocalSupportIds(), entity_set->communicator() ),
            field->dimension() )
    , d_field( field )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Pull data from the application and put it in the vector.
void FieldMultiVector::pullDataFromApplication()
{
    Teuchos::ArrayView<const SupportId> field_supports =
        d_field->getLocalSupportIds();

    int num_supports = field_supports.size();
    int dim = d_field->dimension();

    for ( int d = 0; d < dim; ++d )
    {
        Teuchos::ArrayRCP<double> vector_view = this->getDataNonConst( d );
        for ( int n = 0; n < num_supports; ++n )
        {
            vector_view[n] = d_field->readFieldData( field_supports[n], d );
        }
    }
}

//---------------------------------------------------------------------------//
// Push data from the vector into the application.
void FieldMultiVector::pushDataToApplication()
{
    Teuchos::ArrayView<const SupportId> field_supports =
        d_field->getLocalSupportIds();

    int num_supports = field_supports.size();
    int dim = d_field->dimension();

    for ( int d = 0; d < dim; ++d )
    {
        Teuchos::ArrayRCP<double> vector_view = this->getDataNonConst( d );
        for ( int n = 0; n < num_supports; ++n )
        {
            d_field->writeFieldData( field_supports[n], d, vector_view[n] );
        }
    }

    d_field->finalizeAfterWrite();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_FieldMultiVector.cpp
//---------------------------------------------------------------------------//
