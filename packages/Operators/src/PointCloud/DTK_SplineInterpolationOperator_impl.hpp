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
 * \file   DTK_SplineInterpolationOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Parallel spline interpolator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEINTERPOLATIONOPERATOR_IMPL_HPP
#define DTK_SPLINEINTERPOLATIONOPERATOR_IMPL_HPP

#include "DTK_SplineInterpolationOperator.hpp"
#include "DTK_DBC.hpp"
#include "DTK_CenterDistributor.hpp"
#include "DTK_SplineInterpolationPairing.hpp"
#include "DTK_SplineCoefficientMatrix.hpp"
#include "DTK_SplineEvaluationMatrix.hpp"
#include "DTK_SplineProlongationOperator.hpp"
#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_PredicateComposition.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <BelosPseudoBlockGmresSolMgr.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_DefaultMultipliedLinearOp.hpp>
#include <Thyra_DefaultScaledAdjointLinearOp.hpp>
#include <Thyra_DefaultAddedLinearOp.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Basis,int DIM>
SplineInterpolationOperator<Basis,DIM>::SplineInterpolationOperator(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map,
    const Teuchos::ParameterList& parameters )
    : Base( domain_map, range_map )
    , d_domain_entity_dim( 0 )
    , d_range_entity_dim( 0 )
    , d_use_knn( false )
    , d_knn( 0 )
    , d_radius( 0.0 )
{
    // Determine if we are doing kNN search or radius search.
    if( parameters.isParameter("Type of Search") )
    {
	if ( "Radius" == parameters.get<std::string>("Type of Search") )
	{
	    d_use_knn = false;
	}
	else if ( "Nearest Neighbor" == parameters.get<std::string>("Type of Search") )
	{
	    d_use_knn = true;
	}
	else
	{
	    // Otherwise we got an invalid search type.
	    DTK_INSIST( false );
	}
    }

    // If we are doing kNN support get the number of neighbors.
    if( d_use_knn )
    {
	DTK_REQUIRE( parameters.isParameter("Num Neighbors") );
	d_knn = parameters.get<int>("Num Neighbors");
    }
    
    // Otherwise we are doing the radius search so get the basis radius.
    else
    {
	DTK_REQUIRE( parameters.isParameter("RBF Radius") );
	d_radius = parameters.get<double>("RBF Radius");
    }

    // Get the topological dimension of the domain and range entities. This
    // map will use their centroids for the point cloud.
    if ( parameters.isParameter("Domain Entity Dimension") )
    {
	d_domain_entity_dim = parameters.get<int>("Domain Entity Dimension");
    }
    if ( parameters.isParameter("Range Entity Dimension") )
    {
	d_range_entity_dim = parameters.get<int>("Range Entity Dimension");
    }

    // Get the stratimikos parameters if they exist.
    if ( parameters.isSublist("Stratimikos") )
    {
	d_stratimikos_list = Teuchos::rcp(
	    new Teuchos::ParameterList(parameters.sublist("Stratimikos")) );
    }
}

//---------------------------------------------------------------------------//
// Setup the map operator.
template<class Basis,int DIM>
void SplineInterpolationOperator<Basis,DIM>::setupImpl(
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

    // Prolongation operator.
    Teuchos::RCP<const Root> S;

    // Coefficient matrix polynomial component.
    Teuchos::RCP<const Root> P;

    // Coefficient matrix basis component.
    Teuchos::RCP<const Root> M;

    // Evaluation matrix polynomial component.
    Teuchos::RCP<const Root> Q;

    // Evaluation matrix basis component.
    Teuchos::RCP<const Root> N;

    // Build the concrete operators.
    buildConcreteOperators( domain_space, range_space, S, P, M, Q, N );

    // Create an abstract wrapper for S.
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_range_vector_space_S =
    	Thyra::createVectorSpace<Scalar>( S->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_domain_vector_space_S =
    	Thyra::createVectorSpace<Scalar>( S->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LO,GO> > thyra_S =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,LO,GO>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LO,GO> >(
	thyra_S)->constInitialize( 
	    thyra_range_vector_space_S, thyra_domain_vector_space_S, S );

    // Create an abstract wrapper for P.
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_range_vector_space_P =
    	Thyra::createVectorSpace<Scalar>( P->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_domain_vector_space_P =
    	Thyra::createVectorSpace<Scalar>( P->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LO,GO> > thyra_P =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,LO,GO>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LO,GO> >(
	thyra_P)->constInitialize( 
	    thyra_range_vector_space_P, thyra_domain_vector_space_P, P );

    // Create an abstract wrapper for M.
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_range_vector_space_M =
    	Thyra::createVectorSpace<Scalar>( M->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_domain_vector_space_M =
    	Thyra::createVectorSpace<Scalar>( M->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LO,GO> > thyra_M =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,LO,GO>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LO,GO> >(
	thyra_M)->constInitialize( 
	    thyra_range_vector_space_M, thyra_domain_vector_space_M, M );

    // Create an abstract wrapper for Q.
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_range_vector_space_Q =
    	Thyra::createVectorSpace<Scalar>( Q->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_domain_vector_space_Q =
    	Thyra::createVectorSpace<Scalar>( Q->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LO,GO> > thyra_Q =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,LO,GO>() );
    Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LO,GO> >(
	thyra_Q)->constInitialize( 
	    thyra_range_vector_space_Q, thyra_domain_vector_space_Q, Q );

    // Create an abstract wrapper for N.
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_range_vector_space_N =
    	Thyra::createVectorSpace<Scalar>( N->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_domain_vector_space_N =
    	Thyra::createVectorSpace<Scalar>( N->getDomainMap() );
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LO,GO> > thyra_N =
    	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,LO,GO>() );
        Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LO,GO> >(
	    thyra_N)->constInitialize( 
		thyra_range_vector_space_N, thyra_domain_vector_space_N, N );

    // COUPLING MATRIX ASSEMBLY: A = (Q + N)*[(P + M + P^T)^-1]*S
    // Create a transpose of P.
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyra_P_T =
	Thyra::transpose<Scalar>( thyra_P );

    // Create a composite operator C = (P + M + P^T)
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyra_PpM =
	Thyra::add<Scalar>( thyra_P, thyra_M );
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyra_C =
	Thyra::add<Scalar>( thyra_PpM, thyra_P_T );

    // If we didnt get stratimikos parameters from the input list, create some
    // here.
    if ( Teuchos::is_null(d_stratimikos_list) )
    {
	d_stratimikos_list = Teuchos::parameterList("Stratimikos");

        d_stratimikos_list->set( "Linear Solver Type", "Belos" );
        d_stratimikos_list->set( "Preconditioner Type", "None" );

        auto& linear_solver_types_list =
            d_stratimikos_list->sublist("Linear Solver Types");
        auto& belos_list = linear_solver_types_list.sublist( "Belos" );
        belos_list.set( "Solver Type", "Pseudo Block GMRES" );
        auto& solver_types_list = belos_list.sublist( "Solver Types" );
        auto& gmres_list = solver_types_list.sublist("Pseudo Block GMRES");
        gmres_list.set("Convergence Tolerance", 1.0e-10 );
        gmres_list.set("Verbosity",
                    Belos::Errors + Belos::Warnings +
                    Belos::TimingDetails + Belos::FinalSummary +
                    Belos::StatusTestDetails );
        gmres_list.set("Output Frequency", 1 );
    }

    // Create the inverse of the composite operator C.
    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( d_stratimikos_list );
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory = 
	Thyra::createLinearSolveStrategy( builder );
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyra_C_inv =
    	Thyra::inverse<Scalar>( *factory, thyra_C );

    // Create the composite operator B = (Q + N);
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyra_B =
	Thyra::add<Scalar>( thyra_Q, thyra_N );

    // Create the coupling matrix A = (B * C^-1 * S).
    d_coupling_matrix =
	Thyra::multiply<Scalar>( thyra_B, thyra_C_inv, thyra_S );
    DTK_ENSURE( Teuchos::nonnull(d_coupling_matrix) );
}

//---------------------------------------------------------------------------//
// Apply the operator.
template<class Basis,int DIM>
void SplineInterpolationOperator<Basis,DIM>::applyImpl( 
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
// Transpose apply option.
template<class Basis,int DIM>
bool SplineInterpolationOperator<Basis,DIM>::hasTransposeApplyImpl() const
{
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the concrete operators.
 */
template<class Basis,int DIM>
void SplineInterpolationOperator<Basis,DIM>::buildConcreteOperators(
	const Teuchos::RCP<FunctionSpace>& domain_space,
	const Teuchos::RCP<FunctionSpace>& range_space,
	Teuchos::RCP<const Root>& S,
	Teuchos::RCP<const Root>& P,
	Teuchos::RCP<const Root>& M,
	Teuchos::RCP<const Root>& Q,
	Teuchos::RCP<const Root>& N ) const
{
    // Extract the Support maps.
    const Teuchos::RCP<const typename Base::TpetraMap> domain_map
	= this->getDomainMap();
    const Teuchos::RCP<const typename Base::TpetraMap> range_map
	= this->getRangeMap();

    // Determine if we have range and domain data on this process.
    bool nonnull_domain = Teuchos::nonnull( domain_space->entitySet() );
    bool nonnull_range = Teuchos::nonnull( range_space->entitySet() );

    // Get the parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = domain_map->getComm();

    // Extract the source nodes and their ids.
    Teuchos::ArrayRCP<double> source_centers;
    Teuchos::ArrayRCP<GO> source_support_ids;
    getNodeCoordsAndIds( domain_space, d_domain_entity_dim,
                         source_centers, source_support_ids );
    
    // Extract the target nodes and their ids.    
    Teuchos::ArrayRCP<double> target_centers;
    Teuchos::ArrayRCP<GO> target_support_ids;
    getNodeCoordsAndIds( range_space, d_range_entity_dim,
                         target_centers, target_support_ids );

    // Calculate an approximate neighborhood distance for the local source
    // centers. If using kNN, compute an approximation. If doing a radial
    // search, use the radius.
    double source_proximity = 0.0;
    if ( d_use_knn )
    {
	// Get the local bounding box.
	Teuchos::Tuple<double,6> local_box;
	domain_space->entitySet()->localBoundingBox( local_box );

	// Calculate the largest span of the cardinal directions.
	source_proximity = local_box[3] - local_box[0];
	for ( int d = 1; d < DIM; ++d )
	{
	    source_proximity = std::max( source_proximity,
					 local_box[d+3] - local_box[d] );
	}

	// Take the proximity to be 10% of the largest distance.
	source_proximity *= 0.1;
    }
    else
    {
	source_proximity = d_radius;
    }

    // Gather the source centers that are in the proximity of the source
    // centers on this proc.
    Teuchos::Array<double> dist_sources;
    CenterDistributor<DIM> source_distributor( 
	comm, source_centers(), source_centers(), source_proximity, dist_sources );
    
    // Distribute the global source ids.
    Teuchos::Array<GO> dist_source_support_ids( 
	source_distributor.getNumImports() );
    Teuchos::ArrayView<const GO> source_support_ids_view = source_support_ids();
    Teuchos::ArrayView<GO> dist_gids_view = dist_source_support_ids();
    source_distributor.distribute( source_support_ids_view, dist_gids_view );

    // Build the source/source pairings.
    SplineInterpolationPairing<DIM> source_pairings( 
	dist_sources(), source_centers(), d_use_knn, d_knn, d_radius );

    // Build the basis.
    Teuchos::RCP<Basis> basis = BP::create();

    // PROLONGATION OPERATOR.
    GO offset = comm->getRank() ? 0 : DIM + 1;
    S =	Teuchos::rcp( new SplineProlongationOperator(offset,domain_map) );

    // Get the operator map.
    Teuchos::RCP<const Tpetra::Map<int,GO> > prolongated_map = S->getRangeMap();

    // COEFFICIENT OPERATORS.
    // Build the coefficient operators.
    SplineCoefficientMatrix<Basis,DIM> C( 
	prolongated_map,
	source_centers(), source_support_ids(),
	dist_sources(), dist_source_support_ids(),
	source_pairings, *basis );
    P = C.getP();
    M = C.getM();

    // Cleanup.
    dist_source_support_ids.clear();
    
    // EVALUATION OPERATORS. 
    // Calculate an approximate neighborhood distance for the local target
    // centers. If using kNN, compute an approximation. If doing a radial
    // search, use the radius.
    double target_proximity = 0.0;
    if ( d_use_knn )
    {
	// Get the local bounding box.
	Teuchos::Tuple<double,6> local_box;
	range_space->entitySet()->localBoundingBox( local_box );

	// Calculate the largest span of the cardinal directions.
	target_proximity = local_box[3] - local_box[0];
	for ( int d = 1; d < DIM; ++d )
	{
	    target_proximity = std::max( target_proximity,
					 local_box[d+3] - local_box[d] );
	}

	// Take the proximity to be 10% of the largest distance.
	target_proximity *= 0.1;
    }
    else
    {
	target_proximity = d_radius;
    }

    // Gather the source centers that are in the proximity of the target
    // centers on this proc.
    CenterDistributor<DIM> target_distributor( 
	comm, source_centers(), target_centers(), target_proximity, dist_sources  );

    // Distribute the global source ids.
    dist_source_support_ids.resize( 
	target_distributor.getNumImports() );
    source_support_ids_view = source_support_ids();
    dist_gids_view = dist_source_support_ids();
    target_distributor.distribute( source_support_ids_view, dist_gids_view );

    // Build the source/target pairings.
    SplineInterpolationPairing<DIM> target_pairings( 
	dist_sources(), target_centers(), d_use_knn, d_knn, d_radius );

    // Build the transformation operators.
    SplineEvaluationMatrix<Basis,DIM> B( 
	prolongated_map, range_map,
	target_centers(), target_support_ids(),
	dist_sources(), dist_source_support_ids(),
	target_pairings, *basis );
    N = B.getN();
    Q = B.getQ();
    
    DTK_ENSURE( Teuchos::nonnull(S) );
    DTK_ENSURE( Teuchos::nonnull(P) );
    DTK_ENSURE( Teuchos::nonnull(M) );
    DTK_ENSURE( Teuchos::nonnull(Q) );
    DTK_ENSURE( Teuchos::nonnull(N) );
}

//---------------------------------------------------------------------------//
// Extract node coordinates and ids from an iterator.
template<class Basis,int DIM>
void SplineInterpolationOperator<Basis,DIM>::getNodeCoordsAndIds(
    const Teuchos::RCP<FunctionSpace>& space,    
    const int entity_dim,
    Teuchos::ArrayRCP<double>& centers,
    Teuchos::ArrayRCP<GO>& support_ids ) const
{
    // Get an iterator over the local nodes.
    EntityIterator iterator;
    if ( Teuchos::nonnull(space->entitySet()) )
    {
        LocalEntityPredicate local_predicate(
            space->entitySet()->communicator()->getRank() );	
        PredicateFunction predicate =
            PredicateComposition::And(
                space->selectFunction(),local_predicate.getFunction() );
        iterator = space->entitySet()->entityIterator( entity_dim, predicate );
    }
    
    // Extract the coordinates and support ids of the nodes.
    int local_num_node = iterator.size();
    centers = Teuchos::ArrayRCP<double>( DIM*local_num_node);
    support_ids = Teuchos::ArrayRCP<GO>( local_num_node );
    Teuchos::Array<SupportId> node_supports;
    EntityIterator begin = iterator.begin();
    EntityIterator end = iterator.end();
    int entity_counter = 0;
    for ( EntityIterator entity = begin;
	  entity != end;
	  ++entity, ++entity_counter )
    {
	space->shapeFunction()->entitySupportIds(
	    *entity, node_supports );
	DTK_CHECK( 1 == node_supports.size() );
	support_ids[entity_counter] = node_supports[0];
	space->localMap()->centroid(
	    *entity, centers(DIM*entity_counter,DIM) );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEINTERPOLATIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineInterpolationOperator_impl.hpp
//---------------------------------------------------------------------------//

