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
 * \brief DTK_STKMeshEntity.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHENTITY_HPP
#define DTK_STKMESHENTITY_HPP

#include "DTK_Entity.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Ptr.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshEntity
  \brief STKMesh entity interface definition.
*/
//---------------------------------------------------------------------------//
class STKMeshEntity : public Entity
{
  public:
    /*!
     * \brief Constructor.
     * \param stk_entity The entity to wrap this interface around.
     * \param bulk_data A pointer to the bulk data. We will store a copy of
     * this pointer but not reference count it. We do this because we will
     * create *many* copies of this pointer and do not want to incur the
     * reference counting overhead. We will always make sure that the pointer
     * is in scope both inside and outside of this class while this class
     * exists.
     */
    STKMeshEntity( const stk::mesh::Entity &stk_entity,
                   const Teuchos::Ptr<stk::mesh::BulkData> &bulk_data );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHENTITY_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntity.hpp
//---------------------------------------------------------------------------//
