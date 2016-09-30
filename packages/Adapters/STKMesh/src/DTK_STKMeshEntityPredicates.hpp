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
 * \brief DTK_STKMeshEntityPredicates.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity predicates.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHENTITYPREDICATES_HPP
#define DTK_STKMESHENTITYPREDICATES_HPP

#include <functional>
#include <string>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class STKPartPredicate
 * \part Predicate base class.
 */
class STKPartPredicate
{
  public:

    STKPartPredicate() { /* ... */ }

    bool operator()( Entity entity );

    PredicateFunction getFunction() const
    { return PredicateFunction(*this); }

  protected:

    // Part ids.
    Teuchos::Array<int> b_part_ids;
};

//---------------------------------------------------------------------------//
/*!
  \class STKPartNamePredicate
  \brief Predicates for selecting entities in parts by part name.
*/
class STKPartNamePredicate : public STKPartPredicate
{
  public:

    STKPartNamePredicate( const Teuchos::Array<std::string>& part_names,
                          const Teuchos::RCP<stk::mesh::BulkData>& bulk_data );
};

//---------------------------------------------------------------------------//
/*!
  \class STKPartVectorPredicate
  \brief Predicates for selecting entities in a part vector.
*/
class STKPartVectorPredicate : public STKPartPredicate
{
  public:

    STKPartVectorPredicate( const stk::mesh::PartVector& parts );
};

//---------------------------------------------------------------------------//
/*!
  \class STKSelectorPredicate
  \brief Predicates for selecting entities in a selector.
*/
class STKSelectorPredicate : public STKPartPredicate
{
  public:

    STKSelectorPredicate( const stk::mesh::Selector& selector );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHENTITYPREDICATES_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityPredicates.hpp
//---------------------------------------------------------------------------//
