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
 * \file DTK_TransferTypes.hpp
 * \author Stuart R. Slattery
 * \brief Enumerated types for the transfer subpackage.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TRANSFERTYPES_HPP
#define DTK_TRANSFERTYPES_HPP

namespace DataTransferKit
{

//! Geometry type enumerations.
enum DTK_GeometryType {
    DTK_GeometryType_MIN = 0,
    DTK_COORDINATE_FIELD = DTK_GeometryType_MIN,
    DTK_MESH,
    DTK_GeometryType_MAX = DTK_MESH
};

//! Map type enumerations.
enum DTK_MapType {
    DTK_MapType_MIN = 0,
    DTK_SHARED_DOMAIN = DTK_MapType_MIN,
    DTK_MapType_MAX = DTK_SHARED_DOMAIN
};

} // end namespace DataTransferKit

#endif // end DTK_TRANSFERTYPES_HPP

//---------------------------------------------------------------------------//
// end DTK_TransferTypes.hpp
//---------------------------------------------------------------------------//
