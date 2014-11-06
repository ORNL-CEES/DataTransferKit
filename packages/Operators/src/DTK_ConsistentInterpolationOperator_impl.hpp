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

    // Build the range map for the coupling matrix. 
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_map =
	createDOFMap( range_iterator, range_space->shapeFunction() );

    // Build the domain map for the coupling matrix. 
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_map =
	createDOFMap( domain_iterator, domain_space->shapeFunction() );

    // Allocate the coupling matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,int,std::size_t> > coupling_matrix =
	Tpetra::createCrsMatrix<Scalar,int,std::size_t>( range_map );

    // Construct the entries of the coupling matrix.
    Teuchos::Array<EntityId> domain_entity_ids;
    Teuchos::Array<EntityId>::const_iterator domain_entity_id_it;
    Teuchos::Array<std::size_t> domain_dof_ids;
    Entity domain_entity;
    Teuchos::ArrayView<const double> range_parametric_coords;
    Teuchos::Array<double> domain_shape_values;
    Teuchos::ArrayView<const std::size_t> range_dof_ids = 
	range_map->getNodeElementList();
    Teuchos::ArrayView<const std::size_t>::const_iterator range_dof_id_it;
    EntityIterator range_it;
    EntityIterator range_begin = range_iterator.begin();
    EntityIterator range_end = range_iterator.end();
    for ( range_it = range_begin, range_dof_id_it = range_dof_ids.begin(); 
	  range_it != range_end; 
	  ++range_it, ++range_dof_id_it )
    {
	// Get the domain entities this to which this range entity mapped.
	psearch.getDomainEntitiesFromRange( 
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
	    psearch.rangeParametricCoordinatesInDomain(
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
