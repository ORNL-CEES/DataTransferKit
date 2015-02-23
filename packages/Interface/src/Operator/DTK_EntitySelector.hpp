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
 * \brief DTK_EntitySelector.hpp
 * \author Stuart R. Slattery
 * \brief Entity selector.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYSELECTOR_HPP
#define DTK_ENTITYSELECTOR_HPP

#include <functional>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_EntityIterator.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntitySelector
  \brief Space of a function.

  EntitySelector binds the functional support of a field to a parallel vector
  space.
*/
//---------------------------------------------------------------------------//
class EntitySelector
{
  public:

    /*!
     * \brief Constructor.
     */
    EntitySelector( 
	const EntityType entity_type,
	const std::function<bool(Entity)>& select_function = selectAll );

    /*!
     * \brief Destructor.
     */
    ~EntitySelector();

    /*!
     * \brief Get the entity type to select.
     */
    EntityType entityType() const;

    /*!
     * \brief Get the selector function.
     */
    std::function<bool(Entity)> selectFunction() const;

    /*!
     * \brief Default select function.
     */
    static inline bool selectAll( Entity ) { return true; }

  private:

    // The entity type to select.
    EntityType d_entity_type;

    // The selector function.
    std::function<bool(Entity)> d_select_function;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYSELECTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_EntitySelector.hpp
//---------------------------------------------------------------------------//
