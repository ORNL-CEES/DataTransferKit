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
 * \brief DTK_ConsistentInterpolationOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation oprator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATIONOPERATOR_IMPL_HPP
#define DTK_CONSISTENTINTERPOLATIONOPERATOR_IMPL_HPP

#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "DTK_DBC.hpp"
#include "DTK_ParallelSearch.hpp"
#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_PredicateComposition.hpp"

#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Distributor.hpp>

#include <Thyra_VectorSpaceBase.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraLinearOp.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Scalar>
ConsistentInterpolationOperator<Scalar>::ConsistentInterpolationOperator(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map )
    : Base( domain_map, range_map )
    , d_missed_range_entity_ids( 0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Setup the map operator.
template<class Scalar>
void ConsistentInterpolationOperator<Scalar>::setup(
    const Teuchos::RCP<FunctionSpace>& domain_space,
    const Teuchos::RCP<FunctionSpace>& range_space,
    const Teuchos::RCP<Teuchos::ParameterList>& parameters )
{
    DTK_REQUIRE( Teuchos::nonnull(domain_space) );
    DTK_REQUIRE( Teuchos::nonnull(range_space) );
    DTK_REQUIRE( Teuchos::nonnull(parameters) );

    // Extract the DOF maps.
    const Teuchos::RCP<const typename Base::TpetraMap> domain_map
	= this->getDomainMap();
    const Teuchos::RCP<const typename Base::TpetraMap> range_map
	= this->getRangeMap();

    // Get the parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = domain_map->getComm();

    // Determine if we have range and domain data on this process.
    bool nonnull_domain = Teuchos::nonnull( domain_space->entitySet() );
    bool nonnull_range = Teuchos::nonnull( range_space->entitySet() );

    // Get the physical dimension.
    int physical_dimension = 0;
    if ( nonnull_domain )
    {
	physical_dimension = domain_space->entitySet()->physicalDimension();
    }
    else if ( nonnull_range )
    {
	physical_dimension = range_space->entitySet()->physicalDimension();
    }

    // We will only operate on entities that are locally-owned.
    LocalEntityPredicate local_predicate( comm->getRank() );
    
    // Get an iterator over the domain entities.
    EntityIterator domain_iterator;
    if ( nonnull_domain )
    {
	PredicateFunction domain_predicate =
	    PredicateComposition::And(
		domain_space->entitySelector()->selectFunction(),
		local_predicate.getFunction() );
	domain_iterator = domain_space->entitySet()->entityIterator( 
	    domain_space->entitySelector()->entityType(),
	    domain_predicate );
    }

    // Build a parallel search over the domain.
    ParallelSearch psearch( comm, 
			    physical_dimension,
			    domain_iterator,
			    domain_space->localMap(),
			    *parameters );

    // Get an iterator over the range entities.
    EntityIterator range_iterator;
    if ( nonnull_range )
    {
	PredicateFunction range_predicate =
	    PredicateComposition::And(
		range_space->entitySelector()->selectFunction(),
		local_predicate.getFunction() );
	range_iterator = range_space->entitySet()->entityIterator( 
	    range_space->entitySelector()->entityType(),
	    range_predicate );
    } 

    // Search the domain with the range.
    psearch.search( range_iterator, range_space->localMap(), *parameters );

    // If we are keeping track of range entities that were not mapped, extract
    // them.
    d_missed_range_entity_ids =
	Teuchos::Array<EntityId>( psearch.getMissedRangeEntityIds() );

    // Determine the DOF ids for the range entities found in the local domain
    // on this process and the number of domain entities they were found in
    // globally for averaging.
    std::unordered_map<EntityId,GO> range_dof_id_map;
    std::unordered_map<EntityId,GO> range_dof_count_map;
    {
	// Extract the set of local range entities that were found in domain
	// entities.
	Teuchos::Array<int> export_ranks;
	Teuchos::Array<GO> export_data;
	Teuchos::Array<EntityId> domain_ids;
	Teuchos::Array<EntityId>::const_iterator domain_id_it;
	Teuchos::Array<GO> range_dof_ids;
	EntityIterator range_it;
	EntityIterator range_begin = range_iterator.begin();
	EntityIterator range_end = range_iterator.end();
	for ( range_it = range_begin;
	      range_it != range_end;
	      ++range_it )
	{
	    psearch.getDomainEntitiesFromRange( range_it->id(), domain_ids );

	    for ( domain_id_it = domain_ids.begin();
		  domain_id_it != domain_ids.end();
		  ++domain_id_it )
	    {
		range_space->shapeFunction()->entityDOFIds(
		    *range_it, range_dof_ids );
		DTK_CHECK( 1 == range_dof_ids.size() );

		export_ranks.push_back( 
		    psearch.domainEntityOwnerRank(*domain_id_it) );

		export_data.push_back( range_dof_ids[0] );
		export_data.push_back(
		    Teuchos::as<GO>(range_it->id()) );
		export_data.push_back(
		    Teuchos::as<GO>(domain_ids.size()) );
	    }
	}

	// Communicate the range entity DOF data back to the domain parallel
	// decomposition.
	Tpetra::Distributor range_to_domain_dist( comm );
	int num_import = range_to_domain_dist.createFromSends( export_ranks() );
	Teuchos::Array<GO> import_data( 3*num_import );
	Teuchos::ArrayView<const GO> export_data_view = 
	    export_data();
	range_to_domain_dist.doPostsAndWaits( export_data_view,
					      3,
					      import_data() );

	// Map the range entities to their dof ids.
	for ( int i = 0; i < num_import; ++i )
	{
	    range_dof_id_map.insert(
		std::pair<EntityId,GO>(
		    Teuchos::as<EntityId>(import_data[3*i+1]),
		    import_data[3*i]) );
	    range_dof_count_map.insert(
		std::pair<EntityId,GO>(
		    Teuchos::as<EntityId>(import_data[3*i+1]),
		    import_data[3*i+2]) );
	}
    }

    // Allocate the coupling matrix.
    d_coupling_matrix = 
	Tpetra::createCrsMatrix<Scalar,LO,GO>( range_map );

    // Construct the entries of the coupling matrix.
    Teuchos::Array<EntityId> range_entity_ids;
    Teuchos::Array<EntityId>::const_iterator range_entity_id_it;
    Teuchos::ArrayView<const double> range_parametric_coords;
    Teuchos::Array<double> domain_shape_values;
    Teuchos::Array<double>::iterator domain_shape_it;
    Teuchos::Array<GO> domain_dof_ids;
    EntityIterator domain_it;
    EntityIterator domain_begin = domain_iterator.begin();
    EntityIterator domain_end = domain_iterator.end();
    double scale_val = 0.0;
    for ( domain_it = domain_begin; domain_it != domain_end; ++domain_it )
    {
	// Get the domain DOF ids supporting the domain entity.
	domain_space->shapeFunction()->entityDOFIds( 
	    *domain_it, domain_dof_ids );

	// Get the range entities that mapped into this domain entity.
	psearch.getRangeEntitiesFromDomain( domain_it->id(), range_entity_ids );

	// Sum into the global coupling matrix row for each domain.
	for ( range_entity_id_it = range_entity_ids.begin();
	      range_entity_id_it != range_entity_ids.end();
	      ++range_entity_id_it )
	{
	    // Get the parametric coordinates of the range entity in the
	    // domain entity.
	    psearch.rangeParametricCoordinatesInDomain(
		domain_it->id(), *range_entity_id_it, range_parametric_coords );

	    // Evaluate the shape function at the coordinates.
	    domain_space->shapeFunction()->evaluateValue(
		*domain_it, range_parametric_coords, domain_shape_values );
	    DTK_CHECK( domain_shape_values.size() == domain_dof_ids.size() );

	    // Add the entries as a row in the matrix.
	    if ( range_dof_id_map.count(*range_entity_id_it) )
	    {
		// Compute an average interpolated quantity if the range
		// entity was found in multiple domain entities.
		scale_val = 
		    1.0 / range_dof_count_map.find(*range_entity_id_it)->second;
		for ( domain_shape_it = domain_shape_values.begin();
		      domain_shape_it != domain_shape_values.end();
		      ++domain_shape_it )
		{
		    *domain_shape_it *= scale_val;
		}

		// Consistent interpolation requires one DOF per range
		// entity. Load the row for this range DOF into the matrix.
		d_coupling_matrix->insertGlobalValues( 
		    range_dof_id_map.find(*range_entity_id_it)->second,
		    domain_dof_ids(),
		    domain_shape_values() );
	    }
	}
    }

    // Finalize the coupling matrix.
    d_coupling_matrix->fillComplete( domain_map, range_map );
    DTK_ENSURE( Teuchos::nonnull(d_coupling_matrix) );
}

//---------------------------------------------------------------------------//
// Apply the operator.
template<class Scalar>
void ConsistentInterpolationOperator<Scalar>::applyImpl( 
    const TpetraMultiVector& X,
    TpetraMultiVector &Y,
    Teuchos::ETransp mode,
    Scalar alpha,
    Scalar beta ) const
{
    d_coupling_matrix->apply( X, Y, mode, alpha, beta );
}

//---------------------------------------------------------------------------//
// Return the ids of the range entities that were not mapped during the last
// setup phase (i.e. those that are guaranteed to not receive data from the
// transfer). 
template<class Scalar>
Teuchos::ArrayView<const EntityId> 
ConsistentInterpolationOperator<Scalar>::getMissedRangeEntityIds() const
{
    return d_missed_range_entity_ids();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

# endif // end DTK_CONSISTENTINTERPOLATIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolationOperator_impl.hpp
//---------------------------------------------------------------------------//
