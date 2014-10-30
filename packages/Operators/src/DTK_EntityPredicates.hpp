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
 * \brief DTK_EntityPredicates.hpp
 * \author Stuart R. Slattery
 * \brief Basic entity predicates.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYPREDICATES_HPP
#define DTK_ENTITYPREDICATES_HPP

#include <functional>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityPredicates
  \brief Basic entity predicates.

  A static class of basic entity predicates.
*/
//---------------------------------------------------------------------------//
class EntityPredicates
{
  public:

    /*!
     * \brief Constructor.
     */
    EntityPredicates() { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~EntityPredicates() { /* ... */ }

    //@{
    //! Surface predicate.
    static bool onSurface( Entity entity );
    //@}

    //@{
    //! Block predicate.
    static void setBlockIds( const Teuchos::Array<int>& block_ids );
    static bool inBlocks( Entity entity );
    //@}

    //@{
    //! Boundary predicate.
    static void setBoundaryIds( const Teuchos::Array<int>& boundary_id );
    static bool onBoundaries( Entity entity );
    //@}

    //@{
    //! Owner rank predicate.
    static void setOwnerRank( const int owner_rank );
    static bool hasOwner( Entity entity );
    //@}

    //@{
    //! Entity type predicate.
    static void setEntityType( const EntityType entity_type );
    static bool isEntityType( Entity entity );
    //@}

    private:

    // Current blocks to check.
    static Teuchos::Array<int> d_block_ids;

    // Current boundaries to check.
    static Teuchos::Array<int> d_boundary_ids;

    // Current owner rank to check.
    static int d_owner_rank;

    // Current entity type to check.
    static EntityType d_entity_type;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYPREDICATES_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityPredicates.hpp
//---------------------------------------------------------------------------//
