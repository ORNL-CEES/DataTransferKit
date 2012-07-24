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
 * \file DTK_FieldManager.hpp
 * \author Stuart R. Slattery
 * \brief Field manager declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDMANAGER_HPP
#define DTK_FIELDMANAGER_HPP

#include "DTK_FieldTraits.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class FieldManager
 * \brief Manager object for fields.
 *
 * The field manager manages a field and its parallel decomposition.
 */
//---------------------------------------------------------------------------//
template<class Field>
class FieldManager
{
  public:

    //@{
    //! Typedefs.
    typedef Field                                               field_type;
    typedef FieldTraits<Field>                                  FT;
    typedef Teuchos::Comm<int>                                  CommType;
    typedef Teuchos::RCP<const CommType>                        RCP_Comm;
    //@}

    // Constructor.
    FieldManager( const Field& field, const RCP_Comm& comm );

    // Destructor.
    ~FieldManager();

    //@{
    //! Get the field.
    const Field& field() const
    { return d_field; }

    Field& field()
    { return d_field; }    
    //@}

    //! Get the communicator for the field.
    const RCP_Comm& comm() const
    { return d_comm; }

  private:

    // Field.
    Field d_field;
    
    // Communicator over which the field is defined.
    RCP_Comm d_comm;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_FieldManager_def.hpp"

//---------------------------------------------------------------------------//

#endif // DTK_FIELDMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldManager.hpp
//---------------------------------------------------------------------------//

