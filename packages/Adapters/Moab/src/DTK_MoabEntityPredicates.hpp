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
 * \brief DTK_MoabEntityPredicates.hpp
 * \author Stuart R. Slattery
 * \brief Moab entity predicates.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABMESHENTITYPREDICATES_HPP
#define DTK_MOABMESHENTITYPREDICATES_HPP

#include <functional>
#include <string>
#include <vector>

#include "DTK_Entity.hpp"
#include "DTK_MoabMeshSetIndexer.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_RCP.hpp>

#include <moab/Interface.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class MoabMeshSetPredicate
 * \part Predicate base class.
 */
class MoabMeshSetPredicate
{
  public:
    //! Mesh set predicate. Will return true if a given entity is in the mesh
    //! set.
    MoabMeshSetPredicate( const moab::EntityHandle &mesh_set,
                          const Teuchos::RCP<MoabMeshSetIndexer> &set_indexer );

    //! Mesh set union predicate. Will return true if a given entity is in any
    //! of the given mesh sets. It may be wiser to use moabs set logic to
    //! build a new mesh set of interest and use the basic constructor to
    //! query for set inclusion.
    MoabMeshSetPredicate( const std::vector<moab::EntityHandle> &mesh_sets,
                          const Teuchos::RCP<MoabMeshSetIndexer> &set_indexer );

    bool operator()( Entity entity );

    PredicateFunction getFunction() const { return PredicateFunction( *this ); }

  private:
    // Mesh set ids.
    Teuchos::Array<int> d_set_ids;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABMESHENTITYPREDICATES_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabEntityPredicates.hpp
//---------------------------------------------------------------------------//
