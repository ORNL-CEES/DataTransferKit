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
 * \brief DTK_STKMeshEntityExtraData.hpp
 * \author Stuart R. Slattery
 * \brief Extra data for STK mesh entities.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHENTITYEXTRADATA_HPP
#define DTK_STKMESHENTITYEXTRADATA_HPP

#include "DTK_EntityExtraData.hpp"

#include <stk_mesh/base/Entity.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshEntityExtraData
  \brief A base class for setting extra data with entities.
*/
//---------------------------------------------------------------------------//
class STKMeshEntityExtraData : public EntityExtraData
{
  public:

    STKMeshEntityExtraData( const stk::mesh::Entity stk_entity ) 
	: d_stk_entity( stk_entity )
    { /* ... */ }

    ~STKMeshEntityExtraData() { /* ... */ }

    // STK mesh entity.
    const stk::mesh::Entity d_stk_entity;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_STKMESHENTITYEXTRADATA_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityExtraData.hpp
//---------------------------------------------------------------------------//
