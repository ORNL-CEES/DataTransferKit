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
 * \brief DTK_ParallelSearch.hpp
 * \author Stuart R. Slattery
 * \brief Parallel search.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PARALLELSEARCH_HPP
#define DTK_PARALLELSEARCH_HPP

#include <functional>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_FunctionSpace.hpp"

#include <Teuchos_RCP.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsGraph.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class ParallelSearch
  \brief Parallel search.

  A parallel search finds which domain entities satisfying the domain
  predicate are correctly mapped to the range entities satisfying the range
  predicate. The result of the search is a graph mapping the domain to the
  range and a multivector containing the parametric coordinates of the range
  entity centroids.
*/
//---------------------------------------------------------------------------//
class ParallelSearch
{
  public:

    /*!
     * \brief Constructor.
     */
    ParallelSearch( const Teuchos::RCP<FunctionSpace>& domain_space,
		    const std::function<bool(Entity)>& domain_predicate,
		    const EntityType domain_entity_type,
		    const Teuchos::RCP<FunctionSpace>& range_space,
		    const std::function<bool(Entity)>& range_predicate,
		    const EntityType range_entity_type = ENTITY_TYPE_NODE )

    /*!
     * \brief Destructor.
     */
    ~ParallelSearch();

    /*!
     * \brief Get the domain-to-range entity_graph.
     */
    Teuchos::RCP<Tpetra::CrsGraph<int,EntityId> > graph() const;

    /*!
     * \brief Get the parametric coordinates of the range entities in the
     * domain entities.
     */
    Teuchos::RCP<Tpetra::MultiVector<double,int,EntityId> > 
    rangeCoordinates() const;

  private:

    // Domain-to-range entity graph.
    Teuchos::RCP<Tpetra::CrsGraph<int,EntityId> > d_dtr_graph;

    // Parametric coordinates of the range entities in the domain
    // entities. The map of the vector is the mapped range entities in the
    // decomposition of the domain entities they were found in.
    Teuchos::RCP<Tpetra::MultiVector<double,int,EntityId> > d_parametric_coords;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PARALLELSEARCH_HPP

//---------------------------------------------------------------------------//
// end DTK_ParallelSearch.hpp
//---------------------------------------------------------------------------//
