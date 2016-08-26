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
 * \brief DTK_L2ProjectionOperator.cpp
 * \author Stuart R. Slattery
 * \brief L2 projection operator.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <map>

#include "DTK_L2ProjectionOperator.hpp"
#include "DTK_DBC.hpp"
#include "DTK_ParallelSearch.hpp"
#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_PredicateComposition.hpp"
#include "DTK_IntegrationPoint.hpp"

#include <Teuchos_OrdinalTraits.hpp>

#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Tpetra_Distributor.hpp>

#include <BelosPseudoBlockCGSolMgr.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_DefaultMultipliedLinearOp.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
L2ProjectionOperator::L2ProjectionOperator(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map,
    const Teuchos::ParameterList& parameters )
    : Base( domain_map, range_map )
{
    // Get the integration order.
    const Teuchos::ParameterList& l2_list = parameters.sublist("L2 Projection");
    DTK_REQUIRE( l2_list.isParameter("Integration Order") );
    d_int_order = l2_list.get<int>("Integration Order");
    
    // Get the search list.
    d_search_list = parameters.sublist( "Search" );
}

//---------------------------------------------------------------------------//
// Setup the map operator.
void L2ProjectionOperator::setupImpl(
    const Teuchos::RCP<FunctionSpace>& domain_space,
    const Teuchos::RCP<FunctionSpace>& range_space )
{
    DTK_REQUIRE( Teuchos::nonnull(domain_space) );
    DTK_REQUIRE( Teuchos::nonnull(range_space) );

    // Get the parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	this->getDomainMap()->getComm();

    // Determine if we have range and domain data on this process.
    bool nonnull_domain = Teuchos::nonnull( domain_space->entitySet() );
    bool nonnull_range = Teuchos::nonnull( range_space->entitySet() );

   // Get an iterator over the domain entities.
    EntityIterator domain_iterator;
    if ( nonnull_domain )
    {
	LocalEntityPredicate local_predicate(
	    domain_space->entitySet()->communicator()->getRank() );	
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
	LocalEntityPredicate local_predicate(
	    range_space->entitySet()->communicator()->getRank() );	
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

    // Assemble the coupling matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> > coupling_matrix;
    assembleCouplingMatrix( domain_space, domain_iterator,
			    range_ip_set, coupling_matrix );

    // Create an abstract wrapper for the mass matrix.
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_range_vector_space_M =
    	Thyra::createVectorSpace<double>( mass_matrix->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_domain_vector_space_M =
    	Thyra::createVectorSpace<double>( mass_matrix->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LO,GO> > thyra_M =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,LO,GO>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LO,GO> >(
	thyra_M)->constInitialize( 
	    thyra_range_vector_space_M, thyra_domain_vector_space_M, mass_matrix );

    // Create an abstract wrapper for the coupling matrix.
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_range_vector_space_A =
    	Thyra::createVectorSpace<double>( coupling_matrix->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_domain_vector_space_A =
    	Thyra::createVectorSpace<double>( coupling_matrix->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LO,GO> > thyra_A =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,LO,GO>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LO,GO> >(
	thyra_A)->constInitialize( 
	    thyra_range_vector_space_A, thyra_domain_vector_space_A, coupling_matrix );

    // Create the inverse of the mass matrix. Use the conjugate gradient
    // method to invert the SPD mass matrix.
    Teuchos::RCP<Teuchos::ParameterList> builder_params =
	Teuchos::parameterList("Stratimikos");

    builder_params->set( "Linear Solver Type", "Belos" );
    builder_params->set( "Preconditioner Type", "None" );

    auto& linear_solver_types_list = builder_params->sublist("Linear Solver Types");
    auto& belos_list = linear_solver_types_list.sublist( "Belos" );
    belos_list.set( "Solver Type", "Pseudo Block CG" );
    auto& solver_types_list = belos_list.sublist( "Solver Types" );
    auto& cg_list = solver_types_list.sublist("Pseudo Block CG");
    cg_list.set("Convergence Tolerance", 1.0e-10 );
    cg_list.set("Verbosity",
                Belos::Errors + Belos::Warnings +
                Belos::TimingDetails + Belos::FinalSummary +
                Belos::StatusTestDetails );
    cg_list.set("Output Frequency", 1 );    
    
    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( builder_params );
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > factory = 
	Thyra::createLinearSolveStrategy( builder );
    Teuchos::RCP<const Thyra::LinearOpBase<double> > thyra_M_inv =
    	Thyra::inverse<double>( *factory, thyra_M );

    // Create the projection operator: Op = M^-1 * A.
    d_l2_operator = Thyra::multiply<double>( thyra_M_inv, thyra_A );
    DTK_ENSURE( Teuchos::nonnull(d_l2_operator) );
}

//---------------------------------------------------------------------------//
// Apply the operator.
void L2ProjectionOperator::applyImpl( 
    const TpetraMultiVector& X,
    TpetraMultiVector &Y,
    Teuchos::ETransp mode,
    double alpha,
    double beta ) const
{
    DTK_REQUIRE( Teuchos::NO_TRANS == mode );
    Teuchos::RCP<const Thyra::MultiVectorBase<double> > thyra_X =
	Thyra::createConstMultiVector<double>( Teuchos::rcpFromRef(X) );
    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_Y =
	Thyra::createMultiVector<double>( Teuchos::rcpFromRef(Y) );
    d_l2_operator->apply( 
	Thyra::NOTRANS, *thyra_X, thyra_Y.ptr(), alpha, beta );
}

//---------------------------------------------------------------------------//
// Transpose apply option.
bool L2ProjectionOperator::hasTransposeApplyImpl() const
{
    return false;
}

//---------------------------------------------------------------------------//
// Assemble the mass matrix and range integration point set.
void L2ProjectionOperator::assembleMassMatrix(
    const Teuchos::RCP<FunctionSpace>& range_space,
    EntityIterator range_iterator,
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& mass_matrix,
    Teuchos::RCP<IntegrationPointSet>& range_ip_set )
{
    // Initialize output variables.
    Teuchos::RCP<const Teuchos::Comm<int> > range_comm =
	range_space->entitySet()->communicator();
    mass_matrix = Tpetra::createCrsMatrix<Scalar,LO,GO>( this->getRangeMap() );
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
    range_ip.d_physical_coordinates.resize(
	range_space->entitySet()->physicalDimension() );
    Teuchos::Array<Teuchos::Array<double> > int_points;
    Teuchos::Array<double> int_weights;
    Teuchos::Array<Teuchos::Array<double> > shape_evals;
    Teuchos::Array<SupportId> range_support_ids;
    Teuchos::Array<double> mm_values;
    int range_cardinality = 0;
    int num_ip = 0;
    double temp = 0;

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
	    range_local_map->mapToPhysicalFrame(
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
	    for ( int p = 0; p < num_ip; ++p )
	    {
		temp =
		    range_entity_measure * int_weights[p] * shape_evals[p][ni];
		for ( int nj = 0; nj < range_cardinality; ++nj )
		{
		    mm_values[nj] += temp * shape_evals[p][nj];
		}
	    }
	    mass_matrix->insertGlobalValues( range_support_ids[ni],
					     range_support_ids(),
					     mm_values() );
	}
    }

    // Finalize the mass matrix.
    mass_matrix->fillComplete( this->getRangeMap(), this->getRangeMap() );
    DTK_CHECK( mass_matrix->isFillComplete() );

    // Finalize the integration point set.
    range_ip_set->finalize();
}

//---------------------------------------------------------------------------//
void L2ProjectionOperator::assembleCouplingMatrix(
    const Teuchos::RCP<FunctionSpace>& domain_space,
    EntityIterator domain_iterator,
    const Teuchos::RCP<IntegrationPointSet>& range_ip_set,
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& coupling_matrix )
{
    // Get the parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	domain_space->entitySet()->communicator();

    // Get the physical dimension.
    int physical_dimension = domain_space->entitySet()->physicalDimension();

    // Build a parallel search over the domain.
    ParallelSearch psearch( comm, 
			    physical_dimension,
			    domain_iterator,
			    domain_space->localMap(),
			    d_search_list );

    // Search the domain with the range integration point set.
    EntityIterator ip_iterator = range_ip_set->entityIterator();
    psearch.search( ip_iterator, range_ip_set, d_search_list );

    // Extract the set of local range entities that were found in domain
    // entities.
    int global_max_support = range_ip_set->globalMaxSupportSize();
    int ip_stride = global_max_support + 1;
    Teuchos::Array<int> export_ranks;
    Teuchos::Array<EntityId> export_ip_domain_ids;
    Teuchos::Array<double> export_measures_weights;
    Teuchos::Array<double> export_shape_evals;
    Teuchos::Array<SupportId> export_support_ids;
    Teuchos::Array<EntityId> domain_ids;
    Teuchos::Array<EntityId>::const_iterator domain_id_it;
    EntityIterator ip_it;
    EntityIterator ip_begin = ip_iterator.begin();
    EntityIterator ip_end = ip_iterator.end();
    int num_support = 0;
    for ( ip_it = ip_begin; ip_it != ip_end; ++ip_it )
    {
	// Get the domain entities in which the integration point was found.
	psearch.getDomainEntitiesFromRange( ip_it->id(), domain_ids );

	// Get the current integration point.
	const IntegrationPoint& current_ip =
	    range_ip_set->getPoint( ip_it->id() );
	
	// For each supporting domain entity, pair the integration point id
	// and its support id.
	for ( domain_id_it = domain_ids.begin();
	      domain_id_it != domain_ids.end();
	      ++domain_id_it )
	{
	    // Add the domain rank.
	    export_ranks.push_back( 
		psearch.domainEntityOwnerRank(*domain_id_it) );

	    // Add the ip-domain id pair.
	    export_ip_domain_ids.push_back( ip_it->id() );
	    export_ip_domain_ids.push_back( *domain_id_it );
	    
	    // Add the ip weights times measures scaled by the number of
	    // domains in which ip was found.
	    export_measures_weights.push_back( current_ip.d_owner_measure *
					       current_ip.d_integration_weight /
					       domain_ids.size() );

	    // Add the ip shape function evals and support ids with padding.
	    num_support = current_ip.d_owner_support_ids.size();
	    export_shape_evals.push_back( num_support );
	    export_support_ids.push_back( num_support );
	    for ( int i = 0; i < num_support; ++i )
	    {
		export_shape_evals.push_back(
		    current_ip.d_owner_shape_evals[i] );
		export_support_ids.push_back(
		    current_ip.d_owner_support_ids[i] );
	    }
	    for ( int i = num_support; i < global_max_support; ++i )
	    {
		export_shape_evals.push_back( 0.0 );
		export_support_ids.push_back( 0 );
	    }
	}
    }

    // Communicate the integration points to the domain parallel
    // decomposition.
    Tpetra::Distributor range_to_domain_dist( comm );
    int num_import = range_to_domain_dist.createFromSends( export_ranks() );
    Teuchos::Array<EntityId> import_ip_domain_ids( 2*num_import );
    range_to_domain_dist.doPostsAndWaits(
	export_ip_domain_ids().getConst(), 2, import_ip_domain_ids() );
    Teuchos::Array<double> import_measures_weights( num_import );
    range_to_domain_dist.doPostsAndWaits(
	export_measures_weights().getConst(), 1, import_measures_weights() );
    Teuchos::Array<double> import_shape_evals( ip_stride*num_import );
    range_to_domain_dist.doPostsAndWaits(
	export_shape_evals().getConst(), ip_stride, import_shape_evals() );
    Teuchos::Array<SupportId> import_support_ids( ip_stride*num_import );
    range_to_domain_dist.doPostsAndWaits(
	export_support_ids().getConst(), ip_stride, import_support_ids() );

    // Map the ip-domain pairs to a local id.
    std::map<std::pair<EntityId,EntityId>,int> ip_domain_lid_map;
    std::pair<EntityId,EntityId> lid_pair;
    for ( int n = 0; n < num_import; ++n )
    {
	lid_pair.first = import_ip_domain_ids[2*n];
	lid_pair.second = import_ip_domain_ids[2*n+1];
	ip_domain_lid_map[ lid_pair ] = n;
    }

    // Cleanup before filling the matrix.
    import_ip_domain_ids.clear();
    export_ranks.clear();
    export_ip_domain_ids.clear();
    export_measures_weights.clear();
    export_shape_evals.clear();
    export_support_ids.clear();
    
    // Allocate the coupling matrix.
    coupling_matrix = Tpetra::createCrsMatrix<Scalar,LO,GO>( this->getRangeMap() );

    // Construct the entries of the coupling matrix.
    Teuchos::Array<EntityId> ip_entity_ids;
    Teuchos::Array<EntityId>::const_iterator ip_entity_id_it;
    Teuchos::ArrayView<const double> ip_parametric_coords;
    Teuchos::Array<double> domain_shape_values;
    Teuchos::Array<double> cm_values;
    Teuchos::Array<double>::iterator domain_shape_it;
    Teuchos::Array<GO> domain_support_ids;
    EntityIterator domain_it;
    EntityIterator domain_begin = domain_iterator.begin();
    EntityIterator domain_end = domain_iterator.end();
    int local_id = 0;
    int range_cardinality = 0;
    int domain_cardinality = 0;
    double temp = 0.0;
    int ip_index;
    for ( domain_it = domain_begin; domain_it != domain_end; ++domain_it )
    {
	// Get the domain Support ids supporting the domain entity.
	domain_space->shapeFunction()->entitySupportIds( 
	    *domain_it, domain_support_ids );

	// Get the integration points that mapped into this domain entity.
	psearch.getRangeEntitiesFromDomain( domain_it->id(), ip_entity_ids );

	// Sum into the global coupling matrix row at each integration point.
	for ( ip_entity_id_it = ip_entity_ids.begin();
	      ip_entity_id_it != ip_entity_ids.end();
	      ++ip_entity_id_it )
	{
	    // Get the local data id.
	    lid_pair.first = *ip_entity_id_it;
	    lid_pair.second = domain_it->id();
	    DTK_CHECK( ip_domain_lid_map.count(lid_pair) );
	    local_id = ip_domain_lid_map.find( lid_pair )->second;
	    
	    // Get the parametric coordinates of the integration point in the
	    // domain entity.
	    psearch.rangeParametricCoordinatesInDomain(
		domain_it->id(), *ip_entity_id_it, ip_parametric_coords );

	    // Evaluate the shape function at the integration point
	    // parametric coordinates.
	    domain_space->shapeFunction()->evaluateValue(
		*domain_it, ip_parametric_coords, domain_shape_values );
	    DTK_CHECK( domain_shape_values.size() == domain_support_ids.size() );

	    // Fill the coupling matrix.
	    range_cardinality = import_shape_evals[ip_stride*local_id];
	    domain_cardinality = domain_shape_values.size();
	    cm_values.assign( domain_cardinality, 0.0 );
	    for ( int i = 0; i < range_cardinality; ++i )
	    {
		ip_index = ip_stride*local_id + i + 1;
		temp = import_measures_weights[local_id] *
		       import_shape_evals[ip_index];
		for ( int j = 0; j < domain_cardinality; ++j )
		{
		    cm_values[j] = temp * domain_shape_values[j];
		}
		coupling_matrix->insertGlobalValues(
		    import_support_ids[ip_index],
		    domain_support_ids(),
		    cm_values() );
	    }
	}
    }

    // Finalize the coupling matrix.
    coupling_matrix->fillComplete( this->getDomainMap(), this->getRangeMap() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_L2ProjectionOperator.cpp
//---------------------------------------------------------------------------//
