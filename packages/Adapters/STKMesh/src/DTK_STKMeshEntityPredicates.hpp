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
 * \brief DTK_STKMeshEntityPredicates.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity predicates.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYPREDICATES_HPP
#define DTK_ENTITYPREDICATES_HPP

#include <functional>
#include <string>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <stk_mesh/base/BulkData.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class PartBlockPredicate
  \brief Predicates for selecting entities in a part representing a block
*/
class PartBlockPredicate
{
  public:

    PartBlockPredicate( const Teuchos::Array<std::string>& part_names,
			const Teuchos::RCP<stk::mesh::BulkData>& bulk_data );

    ~PartBlockPredicate() { /* ... */ }

    bool operator()( Entity entity );

    std::function<bool(Entity)> getFunction() const
    { return std::function<bool(Entity)>(*this); }

  private:

    // Part block ids.
    Teuchos::Array<int> d_block_ids;
};

//---------------------------------------------------------------------------//
/*!
  \class PartBoundaryPredicate
  \brief Predicates for selecting entities in a part representing a boundary
*/
class PartBoundaryPredicate
{
  public:

    PartBoundaryPredicate( const Teuchos::Array<std::string>& part_names,
			   const Teuchos::RCP<stk::mesh::BulkData>& bulk_data );

    ~PartBoundaryPredicate() { /* ... */ }

    bool operator()( Entity entity );

    std::function<bool(Entity)> getFunction() const
    { return std::function<bool(Entity)>(*this); }

  private:

    // Part boundary ids.
    Teuchos::Array<int> d_boundary_ids;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYPREDICATES_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityPredicates.hpp
//---------------------------------------------------------------------------//
