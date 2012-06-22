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
 * \file DTK_FieldTools.hpp
 * \author Stuart R. Slattery
 * \brief FieldTools definition
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDTOOLS_HPP
#define DTK_FIELDTOOLS_HPP

#include "DTK_FieldTraits.hpp"
#include <DTK_BoundingBox.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

template<class Field>
class FieldTools
{
  public:

    //@{
    //! Typedefs. 
    typedef Field                           field_type;
    typedef FieldTraits<Field>              FT;
    typedef typename FT::value_type         value_type;
    typedef Teuchos::Comm<int>              CommType;
    typedef Teuchos::RCP<const CommType>    RCP_Comm;
    //@}

    //! Constructor.
    FieldTools()
    { /* ... */ }

    //! Destructor.
    ~FieldTools()
    { /* ... */ }

    // Get the local bounding box for a coordinate field.
    static BoundingBox coordLocalBoundingBox( const Field& field );

    // Get the global bounding box for a coordinate field.
    static BoundingBox coordGlobalBoundingBox( const Field& field,
					       const RCP_Comm& comm );

    // Get the infinity norm of a given field dimension.
    static value_type normInf( const std::size_t dim );

    // Get the L1 norm of a given field dimension.
    static value_type norm1( const std::size_t dim );

    // Get the L2 norm of a given field dimension.
    static value_type norm2( const std::size_t dim );
};

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_FieldTools_def.hpp"

#endif // end DTK_FIELDTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldTools.hpp
//---------------------------------------------------------------------------//

