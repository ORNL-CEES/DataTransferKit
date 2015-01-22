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
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = Box(box_ids[i],comm_rank,i,0.0,0.0,i,1.0,1.0,i+1.0);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[i]);
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

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    Teuchos::Array<Entity> points( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = i + 0.5;
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	points[i] = Point(point_ids[i],comm_rank,point);
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
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*i, point_dofs[i] );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( map_op->getMissedRangeEntityIds().size(), 0 );
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
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = Box(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[i]);
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

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp(
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    Teuchos::Array<Entity> points( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	points[i] = Point(point_ids[i],comm_rank,point);
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
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*point_ids[i], point_dofs[i] );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( map_op->getMissedRangeEntityIds().size(), 0 );
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
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = Box(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[i]);
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

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    Teuchos::Array<Entity> points( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	points[i] = Point(point_ids[i],comm_rank,point);
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
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp(
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val = (comm_rank != comm_size-1) ? 2.0*point_ids[i] : 0.0;
	TEST_EQUALITY( test_val, point_dofs[i] );
    }

    // Check that proc zero had all points not found.
    int num_missed = (comm_rank != comm_size-1) ? 0 : 5;
    Teuchos::Array<EntityId> missed_ids(
	map_op->getMissedRangeEntityIds() );
    TEST_EQUALITY( missed_ids.size(), num_missed );
    std::sort( point_ids.begin(), point_ids.end() );
    std::sort( missed_ids.begin(), missed_ids.end() );
    for ( int i = 0; i < num_missed; ++i )
    {
	TEST_EQUALITY( missed_ids[i], point_ids[i] );
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
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = Box(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[i]);
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

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = ( comm_rank != 0 ) ? 5 : 0;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    Teuchos::Array<Entity> points( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	points[i] = Point(point_ids[i],comm_rank,point);
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
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp(
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*point_ids[i], point_dofs[i] );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( map_op->getMissedRangeEntityIds().size(), 0 );
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
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = Box(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[i]);
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

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 10;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    Teuchos::Array<Entity> points( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = comm_rank*5.0 + i + 0.5;
	points[i] = Point(point_ids[i],comm_rank,point);
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
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val = (5.0*comm_rank+i < num_boxes*comm_size) 
			  ? 2.0*(5.0*comm_rank+i)
			  : 0.0;
	TEST_EQUALITY( test_val, point_dofs[i] );
    }

    // Check that proc zero had some points not found.
    int num_missed = (comm_rank != comm_size-1) ? 0 : 5;
    Teuchos::Array<EntityId> missed_ids(
	map_op->getMissedRangeEntityIds() );
    TEST_EQUALITY( missed_ids.size(), num_missed );
    std::sort( point_ids.begin(), point_ids.end() );
    std::sort( missed_ids.begin(), missed_ids.end() );
    for ( int i = 0; i < num_missed; ++i )
    {
	TEST_EQUALITY( missed_ids[i], point_ids[i+5] );
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
    Teuchos::Array<Entity> boxes( num_boxes );
    box_ids[0] = comm_size - comm_rank - 1;
    box_dofs[0] = 2.0*box_ids[0];
    boxes[0] = Box(box_ids[0],comm_rank,box_ids[0],
		   0.0,0.0,box_ids[0],1.0,1.0,box_ids[0]+1.0);
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[0]);
	
    // Construct a local map for the boxes.
    Teuchos::RCP<EntityLocalMap> domain_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the boxes.
    Teuchos::RCP<EntityShapeFunction> domain_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the boxes.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > domain_dof_map =
	createDOFMap( comm, box_ids() );

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp( 
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 1;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    Teuchos::Array<Entity> points( num_points );
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = comm_rank;
    point_ids[0] = comm_rank;
    point_dofs[0] = 0.0;
    points[0] = Point(point_ids[0],comm_rank,point);
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity( points[0] );

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping. We should get an average on some
    // cores because of the location in multiple domains.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val_1 = (comm_rank != 0) ? 2.0*(comm_rank-1.0) : 0.0;
	double test_val_2 = 2.0*comm_rank;
	TEST_EQUALITY( point_dofs[i], (test_val_1+test_val_2)/2.0 );
    }

    // Check that no missed points were found.
    TEST_EQUALITY( map_op->getMissedRangeEntityIds().size(), 0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, global_missed_range_test )
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
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = Box(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[i]);
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

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp(
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points+1 );
    Teuchos::ArrayRCP<double> point_dofs( num_points+1 );
    Teuchos::Array<Entity> points( num_points+1 );
    for ( int i = 0; i < num_points; ++i )
    {
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = 0.5;
	point[1] = 0.5;
	point[2] = point_ids[i] + 0.5;
	points[i] = Point(point_ids[i],comm_rank,point);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(points[i]);
    }

    // Add a bad point.
    point_ids[5] = num_points*comm_rank + 1000;
    point[0] = -100.0;
    point[1] = 0.0;
    point[2] = 0.0;
    points[5] = Point(point_ids[5],comm_rank,point);
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(points[5]);

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*point_ids[i], point_dofs[i] );
    }

    // Check that the bad point was found.
    Teuchos::ArrayView<const EntityId> missed_range =
	map_op->getMissedRangeEntityIds();
    TEST_EQUALITY( missed_range.size(), 1 );
    TEST_EQUALITY( missed_range[0], 
		   Teuchos::as<EntityId>(num_points*comm_rank + 1000) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ConsistentInterpolationOperator, local_missed_range_test )
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
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = Box(box_ids[i],comm_rank,box_ids[i],
		       box_ids[i],box_ids[i],box_ids[i],
		       box_ids[i]+1.0,box_ids[i]+1.0,box_ids[i]+1.0);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(domain_set)->addEntity(boxes[i]);
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

    // Construct a selector for the boxes.
    Teuchos::RCP<EntitySelector> domain_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_VOLUME) );

    // Construct a function space for the boxes.
    Teuchos::RCP<FunctionSpace> domain_space = Teuchos::rcp(
	new FunctionSpace(domain_set,domain_selector,domain_local_map,domain_shape) );

    // Construct a DOF vector for the boxes.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, boxes(), 1, box_dofs() );

    // RANGE SETUP
    // Make a range entity set.
    Teuchos::RCP<EntitySet> range_set =
	Teuchos::rcp( new BasicEntitySet(comm,3) );
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<std::size_t> point_ids( num_points+1 );
    Teuchos::ArrayRCP<double> point_dofs( num_points+1 );
    Teuchos::Array<Entity> points( num_points+1 );
    for ( int i = 0; i < num_points; ++i )
    {
	point_ids[i] = num_points*comm_rank + i;
	point_dofs[i] = 0.0;
	point[0] = point_ids[i] + 0.5;
	point[1] = point_ids[i] + 0.5;
	point[2] = point_ids[i] + 0.5;
	points[i] = Point(point_ids[i],comm_rank,point);
	Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(points[i]);
    }

    // Add a bad point.
    std::size_t id = num_points*comm_rank;
    point_ids[5] = id + 1000;
    point[0] = id + 0.5;
    point[1] = id + 0.5;
    point[2] = id + 1.5;
    points[5] = Point(point_ids[5],comm_rank,point);
    Teuchos::rcp_dynamic_cast<BasicEntitySet>(range_set)->addEntity(points[5]);

    // Construct a local map for the points.
    Teuchos::RCP<EntityLocalMap> range_local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Construct a shape function for the points.
    Teuchos::RCP<EntityShapeFunction> range_shape =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    // Construct a dof map for the points.
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > range_dof_map =
	createDOFMap( comm, point_ids() );

    // Construct a selector for the points.
    Teuchos::RCP<EntitySelector> range_selector = 
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );

    // Construct a function space for the points.
    Teuchos::RCP<FunctionSpace> range_space = Teuchos::rcp( 
	new FunctionSpace(range_set,range_selector,range_local_map,range_shape) );

    // Construct a DOF vector for the points.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_dofs =
	EntityCenteredDOFVector::pullTpetraMultiVectorFromEntitiesAndView(
	    comm, points(), 1, point_dofs() );

    // MAPPING
    // Create a map.
    Teuchos::RCP<ConsistentInterpolationOperator<double> > map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator<double>(domain_dof_map,range_dof_map) );

    // Setup the map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->set<bool>("Track Missed Range Entities",true);
    map_op->setup( 
	domain_space, range_space, parameters );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

    // Push back to the point dofs.
    EntityCenteredDOFVector::pushTpetraMultiVectorToEntitiesAndView(
	*range_dofs, point_dofs() );

    // Check the results of the mapping.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_EQUALITY( 2.0*point_ids[i], point_dofs[i] );
    }

    // Check that the bad point was found.
    Teuchos::ArrayView<const EntityId> missed_range =
	map_op->getMissedRangeEntityIds();
    TEST_EQUALITY( missed_range.size(), 1 );
    TEST_EQUALITY( missed_range[0], 
		   Teuchos::as<EntityId>(num_points*comm_rank + 1000) );
}

//---------------------------------------------------------------------------//
// end tstConsistentInterpolationOperator.cpp
//---------------------------------------------------------------------------//
