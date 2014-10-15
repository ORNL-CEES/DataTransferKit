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
 * \brief DTK_Entity.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITY_HPP
#define DTK_ENTITY_HPP

#include "DTK_EntityImpl.hpp"
#include "DTK_Types.hpp"
#include "DTK_EntityExtraData.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class Entity
  \brief Geometric entity interface definition.
*/
//---------------------------------------------------------------------------//
class Entity
{
  public:

    /*!
     * \brief Constructor.
     */
    Entity();

    /*!
     * \brief Destructor.
     */
    virtual ~Entity();

    /*!
     * \brief Get the entity type.
     * \return The entity type.
     */
    virtual EntityType entityType() const;

    /*!
     * \brief Get the unique global identifier for the entity.
     * \return A unique global identifier for the entity.
     */
    virtual EntityId id() const;
    
    /*!
     * \brief Get the parallel rank that owns the entity.
     * \return The parallel rank that owns the entity.
     */
    virtual int ownerRank() const;

    /*!
     * \brief Return the physical dimension of the entity.
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    virtual int physicalDimension() const;

    virtual void boundingBox( Teuchos::Tuple<double,6>& bounds ) const;

    /*!
     * \brief Determine if an entity is on the surface of the set.
     */
    virtual bool onSurface() const;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    virtual bool inBlock( const int block_id ) const;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    virtual bool onBoundary( const int boundary_id ) const;

    /*!
     * \brief Get the extra data on the entity.
     */
    virtual Teuchos::RCP<EntityExtraData> extraData() const;

  protected:

    // Geometric entity implementation.
    Teuchos::RCP<EntityImpl> b_entity_impl;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ENTITY_HPP

//---------------------------------------------------------------------------//
// end DTK_Entity.hpp
//---------------------------------------------------------------------------//
