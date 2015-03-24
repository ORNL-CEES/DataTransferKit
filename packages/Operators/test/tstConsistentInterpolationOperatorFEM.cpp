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
 * \file tstConsistentInterpolationOperatorFEM.cpp
 * \author Stuart R. Slattery
 * \brief ConsistentInterpolationOperator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MapOperatorFactory.hpp>
#include <DTK_EntitySelector.hpp>
#include <DTK_FunctionSpace.hpp>
#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_EntityShapeFunction.hpp>
#include <DTK_EntityCenteredShapeFunction.hpp>
#include <DTK_EntityCenteredDOFVector.hpp>
#include <DTK_Box.hpp>
#include <DTK_Point.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>

//---------------------------------------------------------------------------//
// Shape function implementation.
//---------------------------------------------------------------------------//
class TestShapeFunction : public DataTransferKit::EntityShapeFunction
{
  public:

    TestShapeFunction( const int dofs_per_entity ) 
	: d_dofs_per_entity( dofs_per_entity )
    { /* ... */ }

    void entityDOFIds( const DataTransferKit::Entity& entity,
		       Teuchos::Array<std::size_t>& dof_ids ) const
    {
	dof_ids.resize( d_dofs_per_entity );
	for ( int i = 0; i < d_dofs_per_entity; ++i )
	{
	    dof_ids[i] = entity.id()*d_dofs_per_entity + i;
	}
    }

    void evaluateValue( 
	const DataTransferKit::Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<double> & values ) const
    {
	values.assign( d_dofs_per_entity, 1.0 / d_dofs_per_entity );
    }

    void evaluateGradient( 
	const DataTransferKit::Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<Teuchos::Array<double> >& gradients ) const
    { /* ... */ }
	
  private:
    int d_dofs_per_entity;
};

//---------------------------------------------------------------------------//
// DOF map.
//---------------------------------------------------------------------------//
Teuchos::RCP<const Tpetra::Map<int,std::size_t> > createDOFMap(
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::ArrayView<const std::size_t>& entity_ids,
    const int dofs_per_entity )
{
    int num_entity = entity_ids.size();
    Teuchos::Array<std::size_t> dof_ids( num_entity*dofs_per_entity );
    for ( int i = 0; i < num_entity; ++i )
    {
	for ( int n = 0; n < dofs_per_entity; ++n )
	{
	    dof_ids[ i*dofs_per_entity + n ] = entity_ids[i]*dofs_per_entity + n;
	}
    }
    return Tpetra::createNonContigMap<int,std::size_t>( dof_ids, comm );
}

//---------------------------------------------------------------------------//
// DOF vector.
//---------------------------------------------------------------------------//
template<class Scalar>
Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
createTestDOFVector( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::ArrayView<const std::size_t>& entity_ids,
    const std::size_t lda,
    const std::size_t num_vectors )
{
    // Construct a map.
    int num_entity = entity_ids.size();
    int dofs_per_entity = lda / num_entity;
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	createDOFMap( comm, entity_ids, dofs_per_entity );

    // Build a tpetra multivector.
    return Tpetra::createMultiVector<Scalar,int,std::size_t>( map, num_vectors );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, one_to_one_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // DOMAIN SETUP
    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set = 
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = 5;
    int dofs_per_box = 4;
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(box_ids[i],comm_rank,box_ids[i],
		0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0) );
    }

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	createTestDOFVector<double>(
	    comm, box_ids, dofs_per_box*num_boxes, 1 );

    // Make the dofs.
    Teuchos::ArrayRCP<double> box_dofs = domain_dofs->get1dViewNonConst();
    for ( int i = 0; i < num_boxes; ++i )
    {
	for ( int n = 0; n < dofs_per_box; ++n )
	{
	    box_dofs[i*dofs_per_box+n] = 2.0*box_ids[i] + n;
	}
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new TestShapeFunction(dofs_per_box) );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids(), dofs_per_box );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> coords(3);
    Teuchos::Array<Entity> points( num_points );
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	coords[0] = 0.5;
	coords[1] = 0.5;
	coords[2] = point_ids[i] + 0.5;
	points[i] = Point(point_ids[i],comm_rank,coords);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(points[i]);
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids(), 1 );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp(
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points, 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<std::string>("Map Type", "Consistent Interpolation");
    Teuchos::ParameterList& map_list = parameters->sublist("Consistent Interpolation");
    DataTransferKit::MapOperatorFactory<double> factory;
    Teuchos::RCP<DataTransferKit::MapOperator<double> > map_op =
	factory.create( domain_dof_map, range_dof_map, *parameters );

    // Setup the map.
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val = 0.0;
	for ( int n = 0; n < dofs_per_box; ++n )
	{
	    test_val += 2.0*point_ids[i] + n;
	}
	test_val /= dofs_per_box;
	TEST_EQUALITY( test_val, point_dofs[i] );
    }
}

//---------------------------------------------------------------------------//
// end tstConsistentInterpolationOperatorFEM.cpp
//---------------------------------------------------------------------------//
