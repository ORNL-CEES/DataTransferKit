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
 * \brief DTK_ClassicGeometricExtraData.hpp
 * \author Stuart R. Slattery
 * \brief Extra data for basic geometry objects.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSICGEOMETRICENTITYEXTRADATA_HPP
#define DTK_CLASSICGEOMETRICENTITYEXTRADATA_HPP

#include "DTK_EntityExtraData.hpp"

#include <Teuchos_Ptr.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class ClassicGeometricExtraData
  \brief A base class for setting extra data with entities.
*/
//---------------------------------------------------------------------------//
template <class Geometry>
class ClassicGeometricEntityExtraData : public EntityExtraData
{
  public:
    ClassicGeometricEntityExtraData( const Teuchos::Ptr<Geometry> &geometry )
        : d_geometry( geometry )
    { /* ... */
    }

    // Pointer to the classic geometry implementation.
    Teuchos::Ptr<Geometry> d_geometry;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CLASSICGEOMETRICENTITYEXTRADATA_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicGeometricExtraData.hpp
//---------------------------------------------------------------------------//
