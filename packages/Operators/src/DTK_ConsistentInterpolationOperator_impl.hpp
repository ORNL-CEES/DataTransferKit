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
 * \brief DTK_ConsistentInterpolationOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation oprator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATIONOPERATOR_IMPL_HPP
#define DTK_CONSISTENTINTERPOLATIONOPERATOR_IMPL_HPP

#include <unordered_set>
#include <unordered_map>

#include "DTK_DBC.hpp"
#include "DTK_ParallelSearch.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
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
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::RCP<EntitySelector>& domain_selector,
    const Teuchos::RCP<EntitySelector>& range_selector )
    : d_comm( comm )
    , d_domain_selector( domain_selector )
    , d_range_selector( range_selector )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
template<class Scalar>
ConsistentInterpolationOperator<Scalar>::~ConsistentInterpolationOperator()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Setup the map operator.
template<class Scalar>
void ConsistentInterpolationOperator<Scalar>::setup(
    const Teuchos::RCP<FunctionSpace>& domain_space,
    const Teuchos::RCP<FunctionSpace>& range_space,
    const Teuchos::RCP<Teuchos::ParameterList>& parameters )
{
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

    // Get an iterator over the domain entities.
    EntityIterator domain_iterator;
    if ( nonnull_domain )
    {
	domain_iterator = domain_space->entitySet()->entityIterator( 
	    d_domain_selector->entityType(),
	    d_domain_selector->selectFunction() );
    }

    // Build a parallel search over the domain.
    ParallelSearch psearch( d_comm, 
			    physical_dimension,
			    domain_iterator,
			    domain_space->localMap(),
			    *parameters );

    // Get an iterator over the range entities.
    EntityIterator range_iterator;
    if ( nonnull_range )
    {
	range_iterator = range_space->entitySet()->entityIterator( 
	    d_range_selector->entityType(),
	    d_range_selector->selectFunction() );
    } 

    // Build the range map for the coupling matrix. 
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_map =
	createDOFMap( range_iterator, range_space->shapeFunction() );

    // Search the domain with the range.
    psearch.search( range_iterator, range_space->localMap(), *parameters );

    // Determine the DOF ids for the range entities found in the local domain
    // on this process.
    std::unordered_map<EntityId,std::size_t> range_dof_id_map;
    {
	// Extract the set of local range entities that were found in the
	// domain entities on this process.
	std::unordered_set<EntityId> local_range_id_set;
	Teuchos::Array<EntityId> range_entity_ids;
	EntityIterator domain_it;
	EntityIterator domain_begin = domain_iterator.begin();
	EntityIterator domain_end = domain_iterator.end();
	for ( domain_it = domain_begin; domain_it != domain_end; ++domain_it )
	{
	    psearch.getRangeEntitiesFromDomain( 
		domain_it->id(), range_entity_ids );
	    local_range_id_set.insert( 
		range_entity_ids.begin(), range_entity_ids.end() );
	}
	Teuchos::Array<EntityId> local_range_ids( local_range_id_set.size() );
	local_range_ids.assign( 
	    local_range_id_set.begin(), local_range_id_set.end() );
	local_range_id_set.clear();

	// Get the owning ranks of the range entities found in the domain on
	// this process.
	Teuchos::Array<int> local_range_ranks( local_range_ids.size() );
	Teuchos::Array<EntityId>::const_iterator local_range_id_it;
	Teuchos::Array<int>::iterator local_range_rank_it;
	for ( local_range_id_it = local_range_ids.begin(),
	    local_range_rank_it = local_range_ranks.begin();
	      local_range_id_it != local_range_ids.end();
	      ++local_range_id_it, ++local_range_rank_it )
	{
	    *local_range_rank_it = 
		psearch.rangeEntityOwnerRank( *local_range_id_it );
	}

	// Create a communication plan to move the range DOF ids.
	Tpetra::Distributor range_to_domain_dist( d_comm );
	Teuchos::ArrayView<const EntityId> local_range_ids_view =
	    local_range_ids();
	Teuchos::ArrayView<const int> local_range_ranks_view =
	    local_range_ranks();
	Teuchos::Array<EntityId> export_ids;
	Teuchos::Array<int> export_ranks;
	range_to_domain_dist.createFromRecvs( local_range_ids_view,
					      local_range_ranks_view,
					      export_ids,
					      export_ranks );

	// Extract the range dof ids.
	int num_export = export_ids.size();
	Teuchos::Array<std::size_t> export_data( 2*num_export );
	Teuchos::Array<std::size_t> range_dof_ids;
	Entity range_entity;
	for ( int i = 0; i < num_export; ++i )
	{
	    range_space->entitySet()->getEntity( export_ids[i], range_entity );
	    range_space->shapeFunction()->entityDOFIds(
		range_entity, range_dof_ids );
	    DTK_CHECK( 1 == range_dof_ids.size() );
	    export_data[2*i] = range_dof_ids[0];
	    export_data[2*i+1] = Teuchos::as<std::size_t>(export_ids[i]);
	}

	// Redistribute the range entity DOF ids to the domain parallel
	// decomposition.
	int num_import = local_range_ids.size();
	Teuchos::ArrayView<const std::size_t> export_data_view = export_data();
	Teuchos::Array<std::size_t> imported_dof_ids( 2*num_import );
	range_to_domain_dist.doPostsAndWaits( 
	    export_data_view, 2, imported_dof_ids() );

	// Map the range entities to their dof ids.
	for ( int i = 0; i < num_import; ++i )
	{
	    range_dof_id_map.insert(
		std::pair<EntityId,std::size_t>(
		    Teuchos::as<EntityId>(imported_dof_ids[2*i+1]),
		    imported_dof_ids[2*i]) );
	}
    }

    // Allocate the coupling matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,int,std::size_t> > coupling_matrix =
	Tpetra::createCrsMatrix<Scalar,int,std::size_t>( range_map );

    // Construct the entries of the coupling matrix.
    Teuchos::Array<EntityId> range_entity_ids;
    Teuchos::Array<EntityId>::const_iterator range_entity_id_it;
    Teuchos::ArrayView<const double> range_parametric_coords;
    Teuchos::Array<double> domain_shape_values;
    Teuchos::Array<std::size_t> domain_dof_ids;
    EntityIterator domain_it;
    EntityIterator domain_begin = domain_iterator.begin();
    EntityIterator domain_end = domain_iterator.end();
    for ( domain_it = domain_begin; domain_it != domain_end; ++domain_it )
    {
	// Get the domain DOF ids supporting the domain entity.
	domain_space->shapeFunction()->entityDOFIds( 
	    *domain_it, domain_dof_ids );

	// Get the range entities that mapped into this domain entity.
	psearch.getRangeEntitiesFromDomain( domain_it->id(), range_entity_ids );

	// For now, assume a continuous interpolant and use the first domain entity,
	// if there are any, for the shape function evaluation.
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
		// Consistent interpolation requires one DOF per range
		// entity.
		coupling_matrix->insertGlobalValues( 
		    range_dof_id_map.find(*range_entity_id_it)->second,
		    domain_dof_ids(),
		    domain_shape_values() );
	    }
	}
    }

    // Build the domain map for the coupling matrix. 
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_map =
	createDOFMap( domain_iterator, domain_space->shapeFunction() );

    // Finalize the coupling matrix.
    coupling_matrix->fillComplete( domain_map, range_map );

    // Wrap the coupling matrix with the Thyra interface.
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_range_vector_space =
	Thyra::createVectorSpace<Scalar>( coupling_matrix->getRangeMap() );
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyra_domain_vector_space =
	Thyra::createVectorSpace<Scalar>( coupling_matrix->getDomainMap() );
    Teuchos::RCP<Thyra::TpetraLinearOp<Scalar,int,std::size_t> > 
	thyra_coupling_matrix =
	Teuchos::rcp( new Thyra::TpetraLinearOp<Scalar,int,std::size_t>() );
    thyra_coupling_matrix->initialize( thyra_range_vector_space, 
				       thyra_domain_vector_space, 
				       coupling_matrix );

    // Set the coupling matrix with the base class.
    this->b_coupling_matrix = thyra_coupling_matrix;
    DTK_ENSURE( Teuchos::nonnull(this->b_coupling_matrix) );
}

//---------------------------------------------------------------------------//
//! Given an entity iterator and a shape function for those entities,
//! compute the parallel DOF map.
template<class Scalar>
Teuchos::RCP<const Tpetra::Map<int,std::size_t> > 
ConsistentInterpolationOperator<Scalar>::createDOFMap( 
    const EntityIterator& entity_iterator,
    const Teuchos::RCP<EntityShapeFunction>& shape_function ) const
{
    std::unordered_set<std::size_t> entity_dof_set;
    EntityIterator entity_it;
    EntityIterator entity_begin = entity_iterator.begin();
    EntityIterator entity_end = entity_iterator.end();
    Teuchos::Array<std::size_t> entity_dof_ids;
    for ( entity_it = entity_begin; entity_it != entity_end; ++entity_it )
    {
	shape_function->entityDOFIds( *entity_it, entity_dof_ids );
	entity_dof_set.insert( entity_dof_ids.begin(), entity_dof_ids.end() );
    }
    entity_dof_ids.resize( entity_dof_set.size() );
    entity_dof_ids.assign( entity_dof_set.begin(), entity_dof_set.end() );
    entity_dof_set.clear();
    return Tpetra::createNonContigMap<int,std::size_t>( 
	entity_dof_ids(), d_comm );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

# endif // end DTK_CONSISTENTINTERPOLATIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolationOperator_impl.hpp
//---------------------------------------------------------------------------//
