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

#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "DTK_DBC.hpp"
#include "DTK_ParallelSearch.hpp"

#include <Teuchos_OrdinalTraits.hpp>

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

    // Search the domain with the range.
    psearch.search( range_iterator, range_space->localMap(), *parameters );

    // Determine the DOF ids for the range entities found in the local domain
    // on this process and the number of domain entities they were found in
    // globally for averaging. This creates the requirement of uniquely owned
    // domain and range entities input into the map - no ghosts.
    std::unordered_map<EntityId,std::size_t> range_dof_id_map;
    std::unordered_map<EntityId,std::size_t> range_dof_count_map;
    {
	// Extract the set of local range entities that were found in domain
	// entities.
	Teuchos::Array<int> export_ranks;
	Teuchos::Array<std::size_t> export_data;
	Teuchos::Array<EntityId> domain_ids;
	Teuchos::Array<EntityId>::const_iterator domain_id_it;
	Teuchos::Array<std::size_t> range_dof_ids;
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
		    Teuchos::as<std::size_t>(range_it->id()) );
		export_data.push_back(
		    Teuchos::as<std::size_t>(domain_ids.size()) );
	    }
	}

	// Communicate the range entity DOF data back to the domain parallel
	// decomposition.
	Tpetra::Distributor range_to_domain_dist( d_comm );
	int num_import = range_to_domain_dist.createFromSends( export_ranks() );
	Teuchos::Array<std::size_t> import_data( 3*num_import );
	Teuchos::ArrayView<const std::size_t> export_data_view = export_data();
	range_to_domain_dist.doPostsAndWaits( export_data_view,
					      3,
					      import_data() );

	// Map the range entities to their dof ids.
	for ( int i = 0; i < num_import; ++i )
	{
	    range_dof_id_map.insert(
		std::pair<EntityId,std::size_t>(
		    Teuchos::as<EntityId>(import_data[3*i+1]),
		    import_data[3*i]) );
	    range_dof_count_map.insert(
		std::pair<EntityId,std::size_t>(
		    Teuchos::as<EntityId>(import_data[3*i+1]),
		    import_data[3*i+2]) );
	}
    }

    // Allocate the coupling matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,int,std::size_t> > coupling_matrix =
	Tpetra::createCrsMatrix<Scalar,int,std::size_t>( range_space->dofMap() );

    // Construct the entries of the coupling matrix.
    Teuchos::Array<EntityId> range_entity_ids;
    Teuchos::Array<EntityId>::const_iterator range_entity_id_it;
    Teuchos::ArrayView<const double> range_parametric_coords;
    Teuchos::Array<double> domain_shape_values;
    Teuchos::Array<double>::iterator domain_shape_it;
    Teuchos::Array<std::size_t> domain_dof_ids;
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
		DTK_CHECK( range_dof_count_map.count(*range_entity_id_it) );
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
		coupling_matrix->insertGlobalValues( 
		    range_dof_id_map.find(*range_entity_id_it)->second,
		    domain_dof_ids(),
		    domain_shape_values() );
	    }
	}
    }

    // Finalize the coupling matrix.
    coupling_matrix->fillComplete( 
	domain_space->dofMap(), range_space->dofMap() );

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

} // end namespace DataTransferKit

# endif // end DTK_CONSISTENTINTERPOLATIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolationOperator_impl.hpp
//---------------------------------------------------------------------------//
