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

#include "DTK_DBC.hpp"
#include "DTK_ParallelSearch.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraLinearOp.hpp>

#include <Thyra_VectorSpaceBase.hpp>

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
    DTK_REQUIRE( domain_space->entitySet()->physicalDimension() ==
		 range_space->entitySet()->physicalDimension() );

    // Determine if we have range and domain data on this process.
    bool nonnull_domain = Teuchos::nonnull( domain_space->entitySet() );
    bool nonnull_range = Teuchos::nonnull( range_space->entitySet() );

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
			    domain_space->entitySet()->physicalDimension(),
			    domain_iterator,
			    domain_space->localMap(),
			    parameters );

    // Get an iterator over the range entities.
    EntityIterator range_iterator;
    if ( nonnull_range )
    {
	range_iterator = range_space->entitySet()->entityIterator( 
	    d_range_selector->entityType(),
	    d_range_selector->selectFunction() );
    } 

    // Search the domain with the range.
    psearch->search( range_iterator, range_space->localMap(), parameters );

    // Build the row map for the coupling matrix. 
    std::unordered_set<std::size_t> range_dof_set;
    EntityIterator range_it;
    EntityIterator range_begin = range_iterator.begin();
    EntityIterator range_end = range_iterator.end();
    Teuchos::Array<std::size_t> range_dof_ids;
    for ( range_it = range_begin; range_it != range_end; ++range_it )
    {
	// Get the DOFs supporting this range entity.
	range_space->shapeFunction()->entityDOFIds( 
	    *range_it, range_dof_ids );

	// Consistent interpolation requires one DOF per range entity.
	DTK_CHECK( 1 == range_dof_ids.size() );

	// Add to the unique list of DOF ids.
	range_dof_set.insert( range_dof_ids.begin(), range_dof_ids.end() );
    }
    range_dof_ids.resize( range_dof_set.size() );
    std::copy( range_dof_set.begin(), range_dof_set.end(), range_dof_ids );
    range_dof_set.clear();
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > row_map =
	Tpetra::createNonContigMap<int,std::size_t>( range_dof_ids(), d_comm );

    // Allocate the coupling matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,int,std::size_t> > coupling_matrix =
	Tpetra::createCrsMatrix<Scalar,int,std::size_t>( row_map );

    // Construct the entries of the coupling matrix.
    Teuchos::Array<EntityId> domain_entity_ids;
    Teuchos::Array<EntityId>::const_iterator domain_entity_id_it;
    Teuchos::Array<std::size_t> domain_dof_ids;
    Entity domain_entity;
    Teuchos::ArrayView<const double> range_parametric_coords;
    Teuchos::Array<double> domain_shape_values;
    Teuchos::Array<std::size_t>::const_iterator range_dof_id_it;
    for ( range_it = range_begin, range_dof_id_it = range_dof_ids.begin(); 
	  range_it != range_end; 
	  ++range_it, ++range_dof_id_it )
    {
	// Get the domain entities this to which this range entity mapped.
	psearch->getDomainEntitiesFromRange( 
	    range_it->id(), domain_entity_ids );

	// For now, assume a continuous interpolant and use the first entity,
	// if there are any, for the shape function evaluation.
	if ( 0 < domain_entity_ids.size() )
	{
	    // Get the domain entity from its id.
	    domain_space->entitySet()->getEntity( 
		domain_entity_ids[0], domain_entity );

	    // Get the domain DOF ids supporting the domain entity.
	    domain_space->shapeFunction()->entityDOFIds( 
		domain_entity, domain_dof_ids );

	    // Get the parametric coordinates of the range entity in the
	    // domain entity.
	    psearch->rangeParametricCoordinatesInDomain(
		domain_entity.id(), range_it->id(), range_parametric_coords );

	    // Evaluate the shape function at the coordinates.
	    domain_space->shapeFunction()->evaluateValue(
		domain_entity, range_parametric_coords, domain_shape_values );
	    DTK_CHECK( domain_shape_values.size() == domain_dof_ids.size() );

	    // Add the entries as a row in the matrix.
	    coupling_matrix->insertGlobalValues( *range_dof_id_it,
						 domain_dof_ids(),
						 domain_shape_values() );
	}
    }

    // Finalize the coupling matrix.
    coupling_matrix->fillComplete();

    // Wrap the coupling matrix with the Thyra interface.
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyra_space =
	Thyra::createVectorSpace<double>( row_map );
    Teuchos::RCP<Thyra::TpetraLinearOp<double,int,int> > thyra_coupling_matrix =
	Teuchos::rcp( new Thyra::TpetraLinearOp<double,int,int>() );
    thyra_coupling_matrix->initialize( thyra_space, thyra_space, coupling_matrix );

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
