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
 * \brief DTK_L2ProjectionOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief L2 projection operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_L2PROJECTIONOPERATOR_IMPL_HPP
#define DTK_L2PROJECTIONOPERATOR_IMPL_HPP

#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "DTK_DBC.hpp"
#include "DTK_ParallelSearch.hpp"
#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_PredicateComposition.hpp"
#include "DTK_IntegrationPoint.hpp"

#include <Teuchos_OrdinalTraits.hpp>

#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Tpetra_Distributor.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_DefaultMultipliedLinearOp.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Scalar>
L2ProjectionOperator<Scalar>::L2ProjectionOperator(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map,
    const Teuchos::ParameterList& parameters )
    : Base( domain_map, range_map )
{
    // Get the integration order.
    DTK_REQUIRE( parameters.isParameter("Integration Order") );
    d_int_order = parameters.get<int>("Integration Order");
    
    // Get the search list.
    d_search_list = parameters.sublist( "Search" );
}

//---------------------------------------------------------------------------//
// Setup the map operator.
template<class Scalar>
void L2ProjectionOperator<Scalar>::setup(
    const Teuchos::RCP<FunctionSpace>& domain_space,
    const Teuchos::RCP<FunctionSpace>& range_space )
{
    DTK_REQUIRE( Teuchos::nonnull(domain_space) );
    DTK_REQUIRE( Teuchos::nonnull(range_space) );

    // Extract the Support maps.
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
		domain_space->selectFunction(),	local_predicate.getFunction() );
	domain_iterator = domain_space->entitySet()->entityIterator( 
	    domain_space->entitySet()->physicalDimension(),
	    domain_predicate );
    }

    // Get an iterator over the range entities.
    EntityIterator range_iterator;
    if ( nonnull_range )
    {
	PredicateFunction range_predicate =
	    PredicateComposition::And(
		range_space->selectFunction(), local_predicate.getFunction() );
	range_iterator = range_space->entitySet()->entityIterator(
		range_space->entitySet()->physicalDimension(),
		range_predicate );
    } 

    // Assemble the mass matrix over the range entity set.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> > mass_matrix;
    Teuchos::RCP<IntegrationPointSet> range_ip_set;
    assembleMassMatrix( range_space, range_iterator, mass_matrix, range_ip_set );
}

//---------------------------------------------------------------------------//
// Apply the operator.
template<class Scalar>
void L2ProjectionOperator<Scalar>::applyImpl( 
    const TpetraMultiVector& X,
    TpetraMultiVector &Y,
    Teuchos::ETransp mode,
    Scalar alpha,
    Scalar beta ) const
{
    DTK_REQUIRE( Teuchos::NO_TRANS == mode );
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > thyra_X =
	Thyra::createConstMultiVector<Scalar>( Teuchos::rcpFromRef(X) );
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > thyra_Y =
	Thyra::createMultiVector<Scalar>( Teuchos::rcpFromRef(Y) );
    d_coupling_matrix->apply( 
	Thyra::NOTRANS, *thyra_X, thyra_Y.ptr(), alpha, beta );
}

//---------------------------------------------------------------------------//
// Assemble the mass matrix and range integration point set.
template<class Scalar>
void L2ProjectionOperator<Scalar>::buildMassMatrix(
    	const Teuchos::RCP<FunctionSpace>& range_space,
	EntityIterator range_iterator,
	Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& mass_matrix,
	Teuchos::RCP<IntegrationPointSet>& range_ip_set )
{
    // Initialize output variables.
    Teuchos::RCP<const Teuchos::Comm<int> > range_comm =
	range_space->entitySet()->communicator();
    int comm_rank = range_comm->getRank();
    mass_matrix = Tpetra::createCrsMatrix<Scalar,LO,GO>( range_map );
    range_ip_set = Teuchos::rcp( new IntegrationPointSet(range_comm) );

    // Get function space objects.
    Teuchos::RCP<EntityLocalMap> range_local_map = range_space->localMap();
    Teuchos::RCP<EntityShapeFunction> range_shape_function =
	range_space->shapeFunction();
    Teuchos::RCP<EntityIntegrationRule> range_integration_rule =
	range_space->integrationRule();

    // Initialize data for assembly loop.
    double range_entity_measure = 0.0;
    IntegrationPoint range_ip;
    Teuchos::Array<Teuchos::Array<double> > int_points;
    Teuchos::Array<double> int_weights;
    Teuchos::Array<Teuchos::Array<double> > shape_evals;
    Teuchos::Array<SupportId> range_support_ids;
    Teuchos::Array<double> mm_values;
    int range_cardinality = 0;
    int num_ip = 0;

    // Assemble the mass matrix and integration point set.
    EntityIterator range_it;
    EntityIterator range_begin = range_iterator.begin();
    EntityIterator range_end = range_iterator.end();
    for ( range_it = range_begin;
	  range_it != range_end;
	  ++range_it )
    {
	// Get the support ids of the entity.
	range_shape_function->entitySupportIds( *range_it,
						range_support_ids );
	
	// Get the measure of the entity.
	range_entity_measure = range_local_map->measure( *range_it );

	// Get the integration rule.
	range_integration_rule->getIntegrationRule( *range_it,
						    d_int_order,
						    int_points,
						    int_weights );

	// Create new integration points.
	num_ip = int_weights.size();
	shape_evals.resize( num_ip );
	for ( int p = 0; p < num_ip; ++p )
	{
	    // Add owner data.
	    range_ip.d_owner_gid = range_it->id();
	    range_ip.d_owner_rank = comm_rank;
	    range_ip.d_owner_measure = range_entity_measure;
	    range_ip.d_owner_support_ids = range_support_ids;
	    range_ip.d_integration_weight = int_weights[p];

	    // Evaluate the shape function.
	    range_shape_function->evaluateValue( *range_it,
						 int_points[p](),
						 shape_evals[p] );
	    range_ip.d_owner_shape_evals = shape_evals[p];

	    // Map the integration point to the physical frame of the range
	    // entity.
	    range_ip.d_physical_coordinates.resize(
		range_it->physicalDimension() );
	    range_local_map.mapToPhysicalFrame(
		*range_it,
		int_points[p](),
		range_ip.d_physical_coordinates() );

	    // Add the integration point to the set.
	    range_ip_set->addPoint( range_ip );
	}

	// Fill the mass matrix.
	range_cardinality = range_support_ids.size();
	for ( int ni = 0; ni < range_cardinality; ++ni )
	{
	    mm_values.assign( range_cardinality, 0.0 );
	    for ( int nj = 0; nj < range_cardinality; ++nj )
	    {
		for ( int p = 0; p < num_ip; ++p )
		{
		    mm_values[nj] += range_entity_measure *
				     int_weights[p] *
				     shape_evals[p][ni] *
				     shape_evals[p][nj];
		}
	    }
	    mass_matrix->insertGlobalValues( range_support_ids[ni],
					     range_support_ids(),
					     mm_values() );
	}
    }

    // Finalize the mass matrix.
    mass_matrix->fillComplete();
    DTK_CHECK( mass_matrix->isFillComplete() );

    // Finalize the integration point set.
    range_ip_set->finalize();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

# endif // end DTK_L2PROJECTIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_L2ProjectionOperator_impl.hpp
//---------------------------------------------------------------------------//
