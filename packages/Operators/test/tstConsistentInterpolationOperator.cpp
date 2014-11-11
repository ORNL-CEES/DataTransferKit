//---------------------------------------------------------------------------//
/*! 
 * \file tstConsistentInterpolationOperator.cpp
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

#include <DTK_ConsistentInterpolationOperator.hpp>
#include <DTK_EntitySelector.hpp>
#include <DTK_FunctionSpace.hpp>
#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
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

//---------------------------------------------------------------------------//
// DOF map.
//---------------------------------------------------------------------------//
Teuchos::RCP<const Tpetra::Map<int,std::size_t> > createDOFMap(
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::ArrayView<const std::size_t>& entity_ids )
{
    return Tpetra::createNonContigMap<int,std::size_t>( entity_ids, comm );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, all_to_one_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();

    // DOMAIN SETUP
    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set = 
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_boxes = (comm->getRank() == 0) ? 5 : 0;
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = i;
	box_dofs[i] = 2.0*box_ids[i];
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(box_ids[i],comm_rank,i,0.0,0.0,i,1.0,1.0,i+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids() );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_local_map,domain_shape,domain_dof_map) );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > domain_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, box_ids, box_dofs, num_boxes, 1 );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = i + 0.5;
	bool on_surface = (i%2==0);
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(point_ids[i],comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_local_map,range_shape,range_dof_map) );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > range_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, point_ids, point_dofs, num_points, 1 );

    // MAPPING
    // Create a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(
	    comm,domain_selector,range_selector) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( Thyra::NOTRANS, *domain_dofs, range_dofs.ptr(), 1.0, 0.0 );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*i, point_dofs[i] );
    }
}

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
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(box_ids[i],comm_rank,box_ids[i],
		0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids() );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp(
	new FunctionSpace(domain_set,domain_local_map,domain_shape,domain_dof_map) );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > domain_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, box_ids, box_dofs, num_boxes, 1 );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	bool on_surface = (i%2==0);
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(point_ids[i],comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_local_map,range_shape,range_dof_map) );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > range_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, point_ids, point_dofs, num_points, 1 );

    // MAPPING
    // Create a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(
	    comm,domain_selector,range_selector) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( Thyra::NOTRANS, *domain_dofs, range_dofs.ptr(), 1.0, 0.0 );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*point_ids[i], point_dofs[i] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, no_domain_0_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();
    int comm_rank = comm->getRank();

    // DOMAIN SETUP
    // Make a domain entity set.
    Teuchos::RCP<EntitySet> domain_set = 
	Teuchos::rcp( new BasicEntitySet(comm,3) );

    // Don't put domain entities on proc 0.
    int num_boxes = (comm->getRank() != 0) ? 5 : 0;
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(box_ids[i],comm_rank,box_ids[i],
		0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids() );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_local_map,domain_shape,domain_dof_map) );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > domain_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, box_ids, box_dofs, num_boxes, 1 );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	bool on_surface = (i%2==0);
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(point_ids[i],comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp(
	new FunctionSpace(range_set,range_local_map,range_shape,range_dof_map) );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > range_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, point_ids, point_dofs, num_points, 1 );

    // MAPPING
    // Create a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(
	    comm,domain_selector,range_selector) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( Thyra::NOTRANS, *domain_dofs, range_dofs.ptr(), 1.0, 0.0 );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val = (comm_rank != comm_size-1) ? 2.0*point_ids[i] : 0.0;
	TEST_EQUALITY( test_val, point_dofs[i] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, no_range_0_test )
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
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(box_ids[i],comm_rank,box_ids[i],
		0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids() );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_local_map,domain_shape,domain_dof_map) );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > domain_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, box_ids, box_dofs, num_boxes, 1 );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = ( comm_rank != 0 ) ? 5 : 0;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	bool on_surface = (i%2==0);
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(point_ids[i],comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp(
	new FunctionSpace(range_set,range_local_map,range_shape,range_dof_map) );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > range_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, point_ids, point_dofs, num_points, 1 );

    // MAPPING
    // Create a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(
	    comm,domain_selector,range_selector) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( Thyra::NOTRANS, *domain_dofs, range_dofs.ptr(), 1.0, 0.0 );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*point_ids[i], point_dofs[i] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, many_to_many_test )
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
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	    Box(box_ids[i],comm_rank,box_ids[i],
		0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0) );
    }

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids() );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_local_map,domain_shape,domain_dof_map) );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > domain_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, box_ids, box_dofs, num_boxes, 1 );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 10;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	bool on_surface = (i%2==0);
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = comm_rank*5.0 + i + 0.5;
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	    Point(point_ids[i],comm_rank,point,on_surface) );
    }

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_local_map,range_shape,range_dof_map) );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > range_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, point_ids, point_dofs, num_points, 1 );

    // MAPPING
    // Create a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(
	    comm,domain_selector,range_selector) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( Thyra::NOTRANS, *domain_dofs, range_dofs.ptr(), 1.0, 0.0 );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val = (5.0*comm_rank+i < num_boxes*comm_size) 
			  ? 2.0*(5.0*comm_rank+i)
			  : 0.0;
	TEST_EQUALITY( test_val, point_dofs[i] );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, point_multiple_neighbors_test )
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
    int num_boxes = 1;
    Teuchos::Array<std::size_t> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    box_ids[0] = comm_size - comm_rank - 1;
    box_dofs[0] = 2.0*box_ids[0];
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(
	Box(box_ids[0],box_ids[0],box_ids[0],
	    0.0,0.0,box_ids[0],1.0,1.0,box_ids[0]+1.0) );

    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids() );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_local_map,domain_shape,domain_dof_map) );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > domain_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, box_ids, box_dofs, num_boxes, 1 );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 1;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = comm_rank;
    point_ids[0] = comm_rank;
    point_dofs[0] = 0.0;
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(
	Point(point_ids[0],comm_rank,point,false) );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_local_map,range_shape,range_dof_map) );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > range_dofs =
	EntityCenteredDOFVector::createThyraMultiVector(
	    comm, point_ids, point_dofs, num_points, 1 );

    // MAPPING
    // Create a map.
    Teuchos::RCP<MapOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(
	    comm,domain_selector,range_selector) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    map_op->setup( domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( Thyra::NOTRANS, *domain_dofs, range_dofs.ptr(), 1.0, 0.0 );

    // Check the results of the mapping. We should get an average on some
    // cores because of the location in multiple domains.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val_1 = (comm_rank != 0) ? 2.0*(comm_rank-1.0) : 0.0;
	double test_val_2 = 2.0*comm_rank;
	TEST_EQUALITY( point_dofs[i], (test_val_1+test_val_2)/2.0 );
    }
}

//---------------------------------------------------------------------------//
// end tstConsistentInterpolationOperator.cpp
//---------------------------------------------------------------------------//
