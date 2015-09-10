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
 * \brief DTK_ParallelSearch.hpp
 * \author Stuart R. Slattery
 * \brief Parallel search.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PARALLELSEARCH_HPP
#define DTK_PARALLELSEARCH_HPP

#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_EntityIterator.hpp"
#include "DTK_EntityLocalMap.hpp"
#include "DTK_CoarseGlobalSearch.hpp"
#include "DTK_CoarseLocalSearch.hpp"
#include "DTK_FineLocalSearch.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class ParallelSearch
  \brief Parallel search.

  A parallel search finds which domain entities satisfying the domain
  predicate are correctly mapped to the range entities satisfying the range
  predicate. The result of the search is a graph mapping the domain to the
  range and the parametric coordinates of the range entity centroids in their
  respective domain entities.

  The search has two simultaneous states: one in the parallel decomposition of
  the domain and one in the parallel decomposition of the range. The interface
  functions assume one decomposition or the other.
*/
//---------------------------------------------------------------------------//
class ParallelSearch
{
  public:

    /*!
     * \brief Constructor.
     */
    ParallelSearch( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		    const int physical_dimension,
		    const EntityIterator& domain_iterator,
		    const Teuchos::RCP<EntityLocalMap>& domain_local_map,
		    const Teuchos::ParameterList& parameters );

    /*
     * \brief Search the domain with the range entity centroids and construct
     * the graph. This will update the state of the object.
     */
    void search( const EntityIterator& range_iterator,
		 const Teuchos::RCP<EntityLocalMap>& range_local_map,
		 const Teuchos::ParameterList& parameters );

    /*!
     * \brief Given a domain entity id on a domain process, get the ids of the
     * range entities that mapped to it.
     */
    void getRangeEntitiesFromDomain( 
	const EntityId domain_id, Teuchos::Array<EntityId>& range_ids ) const;

    /*!
     * \brief Given a range entity id on a range process, get the ids of the
     * domain entities that it mapped to.
     */
    void getDomainEntitiesFromRange( 
	const EntityId range_id, Teuchos::Array<EntityId>& domain_ids ) const;

    /*!
     * \brief Get the owner rank of a given range entity on a domain process.
     */
    int rangeEntityOwnerRank( const EntityId range_id ) const;

    /*!
     * \brief Get the owner rank of a given domain entity on a rank process.
     */
    int domainEntityOwnerRank( const EntityId domain_id ) const;

    /*!
     * \brief Get the parametric coordinates of the range entities in the
     * domain entities on a domain process.
     */
    void rangeParametricCoordinatesInDomain( 
	const EntityId domain_id,
	const EntityId range_id,
	Teuchos::ArrayView<const double>& parametric_coords ) const;

    /*!
     * \brief Return the ids of the range entities that were not during the
     * last search (i.e. those that are guaranteed to not receive data from
     * the transfer) on a range process.
     *
     * \return A view of the ids.
     */
    Teuchos::ArrayView<const EntityId> getMissedRangeEntityIds() const;

  private:

    // Parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Phyiscal dimension.
    int d_physical_dim;

    // Empty domain flag.
    bool d_empty_domain;

    // Empty range flag.
    bool d_empty_range;

    // Coarse global search.
    Teuchos::RCP<CoarseGlobalSearch> d_coarse_global_search;

    // Coarse local search.
    Teuchos::RCP<CoarseLocalSearch> d_coarse_local_search;

    // Fine local search.
    Teuchos::RCP<FineLocalSearch> d_fine_local_search;

    // Range owner rank map.
    std::unordered_map<EntityId,int> d_range_owner_ranks;

    // Domain owner rank map.
    std::unordered_map<EntityId,int> d_domain_owner_ranks;

    // Domain-to-range entity map.
    std::unordered_multimap<EntityId,EntityId> d_domain_to_range_map;

    // Range-to-domain entity map.
    std::unordered_multimap<EntityId,EntityId> d_range_to_domain_map;

    // Parametric coordinates of the range entities in the domain
    // entities. The first key is the range id, the second key is the domain
    // id.
    std::unordered_map<EntityId,
		       std::unordered_map<EntityId,Teuchos::Array<double> > 
		       > d_parametric_coords;

    // Boolean for tracking missed range entities.
    bool d_track_missed_range_entities;

    // An array of range entity ids that were not mapped during the last call
    // to setup.
    mutable Teuchos::Array<EntityId> d_missed_range_entity_ids;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PARALLELSEARCH_HPP

//---------------------------------------------------------------------------//
// end DTK_ParallelSearch.hpp
//---------------------------------------------------------------------------//
