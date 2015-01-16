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
 * \file   DTK_SplineInterpolationOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Parallel spline interpolator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEINTERPOLATIONOPERATOR_IMPL_HPP
#define DTK_SPLINEINTERPOLATIONOPERATOR_IMPL_HPP

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

#include <Tpetra_Map.hpp>

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
template<class Scalar,class Basis,int DIM>
SplineInterpolationOperator<Scalar,Basis,DIM>::SplineInterpolationOperator(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map )
    : Base( domain_map, range_map )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
template<class Scalar,class Basis,int DIM>
SplineInterpolationOperator<Scalar,Basis,DIM>::~SplineInterpolationOperator()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Setup the map operator.
template<class Scalar,class Basis,int DIM>
void SplineInterpolationOperator<Scalar,Basis,DIM>::setup(
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

    // Make sure we are applying the map to nodes.
    DTK_REQUIRE( domain_space->entitySelector()->entityType() ==
		 ENTITY_TYPE_NODE );
    DTK_REQUIRE( range_space->entitySelector()->entityType() ==
		 ENTITY_TYPE_NODE );

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
    buildConcreteOperators( domain_space, range_space, parameters, S, P, M, Q, N );

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

    // Create the inverse of the composite operator C.
    Teuchos::RCP<Teuchos::ParameterList> builder_params =
      sublist(parameters, "Stratimikos");
    Teuchos::updateParametersFromXmlString(
      "<ParameterList name=\"Stratimikos\">"
        "<Parameter name=\"Linear Solver Type\" type=\"string\" value=\"Belos\"/>"
        "<Parameter name=\"Preconditioner Type\" type=\"string\" value=\"None\"/>"
      "</ParameterList>"
      ,
      builder_params.ptr()
      );
    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( builder_params );
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
template<class Scalar,class Basis,int DIM>
void SplineInterpolationOperator<Scalar,Basis,DIM>::apply( 
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
/*!
 * \brief Build the concrete operators.
 */
template<class Scalar,class Basis,int DIM>
void SplineInterpolationOperator<Scalar,Basis,DIM>::buildConcreteOperators(
	const Teuchos::RCP<FunctionSpace>& domain_space,
	const Teuchos::RCP<FunctionSpace>& range_space,
	const Teuchos::RCP<Teuchos::ParameterList>& parameters,
	Teuchos::RCP<const Root>& S,
	Teuchos::RCP<const Root>& P,
	Teuchos::RCP<const Root>& M,
	Teuchos::RCP<const Root>& Q,
	Teuchos::RCP<const Root>& N ) const
{
    // Get the basis radius.
    DTK_REQUIRE( parameters->isParameter("RBF Radius") );
    double radius = parameters->get<double>("RBF Radius");

    // Extract the DOF maps.
    const Teuchos::RCP<const typename Base::TpetraMap> domain_map
	= this->getDomainMap();
    const Teuchos::RCP<const typename Base::TpetraMap> range_map
	= this->getRangeMap();

    // Determine if we have range and domain data on this process.
    bool nonnull_domain = Teuchos::nonnull( domain_space->entitySet() );
    bool nonnull_range = Teuchos::nonnull( range_space->entitySet() );

    // Get the parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = domain_map->getComm();

    // We will only operate on entities that are locally-owned.
    LocalEntityPredicate local_predicate( comm->getRank() );

    // Extract the source centers and their ids.
    EntityIterator domain_iterator;
    if ( nonnull_domain )
    {
	std::function<bool(Entity)> domain_predicate =
	    PredicateComposition::And(
		domain_space->entitySelector()->selectFunction(),
		local_predicate.getFunction() );
	domain_iterator = domain_space->entitySet()->entityIterator( 
	    domain_space->entitySelector()->entityType(),
	    domain_predicate );
    }
    int local_num_src = domain_iterator.size();
    Teuchos::ArrayRCP<double> source_centers( DIM*local_num_src);
    Teuchos::ArrayRCP<GO> source_dof_ids( local_num_src );
    Teuchos::Array<std::size_t> source_node_dofs;
    EntityIterator domain_begin = domain_iterator.begin();
    EntityIterator domain_end = domain_iterator.end();
    int entity_counter = 0;
    for ( EntityIterator domain_entity = domain_begin;
	  domain_entity != domain_end;
	  ++domain_entity, ++entity_counter )
    {
	domain_space->shapeFunction()->entityDOFIds( *domain_entity, source_node_dofs );
	DTK_CHECK( 1 == source_node_dofs.size() );
	source_dof_ids[entity_counter] = source_node_dofs[0];
	domain_space->localMap()->centroid(
	    *domain_entity, source_centers(DIM*entity_counter,DIM) );
    }

    // Extract the target centers and their ids.
    EntityIterator range_iterator;
    if ( nonnull_range )
    {
	std::function<bool(Entity)> range_predicate =
	    PredicateComposition::And(
		range_space->entitySelector()->selectFunction(),
		local_predicate.getFunction() );
	range_iterator = range_space->entitySet()->entityIterator( 
	    range_space->entitySelector()->entityType(),
	    range_predicate );
    } 
    int local_num_tgt = range_iterator.size();
    Teuchos::ArrayRCP<double> target_centers( DIM*local_num_tgt );
    Teuchos::ArrayRCP<GO> target_dof_ids( local_num_tgt );
    Teuchos::Array<std::size_t> target_node_dofs;
    EntityIterator range_begin = range_iterator.begin();
    EntityIterator range_end = range_iterator.end();
    entity_counter = 0;
    for ( EntityIterator range_entity = range_begin;
	  range_entity != range_end;
	  ++range_entity, ++entity_counter )
    {
	range_space->shapeFunction()->entityDOFIds( *range_entity, target_node_dofs );
	DTK_CHECK( 1 == target_node_dofs.size() );
	target_dof_ids[entity_counter] = target_node_dofs[0];
	range_space->localMap()->centroid(
	    *range_entity, target_centers(DIM*entity_counter,DIM) );
    }

    // PROLONGATION OPERATOR.
    GO offset = comm->getRank() ? 0 : DIM + 1;
    S =	Teuchos::rcp( 
	new SplineProlongationOperator<Scalar,GO>(offset,domain_map) );

    // COEFFICIENT OPERATORS.
    // Gather the source centers that are within a radius of the source
    // centers on this proc.
    Teuchos::Array<double> dist_sources;
    CenterDistributor<DIM> source_distributor( 
	comm, source_centers(), source_centers(), radius, dist_sources );
    
    // Distribute the global source ids.
    Teuchos::Array<GO> dist_source_dof_ids( 
	source_distributor.getNumImports() );
    Teuchos::ArrayView<const GO> source_dof_ids_view = source_dof_ids();
    Teuchos::ArrayView<GO> dist_gids_view = dist_source_dof_ids();
    source_distributor.distribute( source_dof_ids_view, dist_gids_view );

    // Build the source/source pairings.
    SplineInterpolationPairing<DIM> source_pairings( 
	dist_sources(), source_centers(), radius );

    // Build the basis.
    Teuchos::RCP<Basis> basis = BP::create( radius );

    // Get the operator map.
    Teuchos::RCP<const Tpetra::Map<int,GO> > prolongated_map = S->getRangeMap();

    // Build the coefficient operators.
    SplineCoefficientMatrix<Basis,DIM> C( 
	prolongated_map,
	source_centers(), source_dof_ids(),
	dist_sources(), dist_source_dof_ids(),
	source_pairings, *basis );
    P = C.getP();
    M = C.getM();

    // Cleanup.
    dist_source_dof_ids.clear();
    
    // EVALUATION OPERATORS. 
    // Gather the source centers that are within a radius of the target
    // centers on this proc.
    CenterDistributor<DIM> target_distributor( 
	comm, source_centers(), target_centers(), radius, dist_sources  );

    // Distribute the global source ids.
    dist_source_dof_ids.resize( 
	target_distributor.getNumImports() );
    source_dof_ids_view = source_dof_ids();
    dist_gids_view = dist_source_dof_ids();
    target_distributor.distribute( source_dof_ids_view, dist_gids_view );

    // Build the source/target pairings.
    SplineInterpolationPairing<DIM> target_pairings( 
	dist_sources(), target_centers(), radius );

    // Build the transformation operators.
    SplineEvaluationMatrix<Basis,DIM> B( 
	prolongated_map, range_map,
	target_centers(), target_dof_ids(),
	dist_sources(), dist_source_dof_ids(),
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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEINTERPOLATIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineInterpolationOperator_impl.hpp
//---------------------------------------------------------------------------//

