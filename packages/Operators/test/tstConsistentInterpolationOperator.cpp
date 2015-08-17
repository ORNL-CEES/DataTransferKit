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
#include <DTK_BasicGeometryManager.hpp>
#include <DTK_EntityCenteredField.hpp>
#include <DTK_FieldMultiVector.hpp>
#include <DTK_BoxGeometry.hpp>
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
    int num_boxes = (comm->getRank() == 0) ? 5 : 0;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = BoxGeometry(box_ids[i],comm_rank,i,0.0,0.0,i,1.0,1.0,i+1.0);
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // RANGE SETUP
    // Make a range entity set.
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
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
    }

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup(
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
    int num_boxes = 5;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = BoxGeometry(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // RANGE SETUP
    // Make a range entity set.
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
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
    }

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup(
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
    // Don't put domain entities on proc 0.
    int num_boxes = (comm->getRank() != 0) ? 5 : 0;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = BoxGeometry(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );    

    // RANGE SETUP
    // Make a range entity set.
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
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
    }

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup( 
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
    int num_boxes = 5;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = BoxGeometry(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // RANGE SETUP
    // Make a range entity set.
    int num_points = ( comm_rank != 0 ) ? 5 : 0;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
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
    }

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);    
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup( 
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
    int num_boxes = 5;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = BoxGeometry(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // RANGE SETUP
    // Make a range entity set.
    int num_points = 10;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
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
    }

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup( 
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
    int num_boxes = 1;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    box_ids[0] = comm_size - comm_rank - 1;
    box_dofs[0] = 2.0*box_ids[0];
    boxes[0] = BoxGeometry(box_ids[0],comm_rank,box_ids[0],
		   0.0,0.0,box_ids[0],1.0,1.0,box_ids[0]+1.0);
	
    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // RANGE SETUP
    // Make a range entity set.
    int num_points = 1;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points );
    Teuchos::ArrayRCP<double> point_dofs( num_points );
    Teuchos::Array<Entity> points( num_points );
    point[0] = 0.5;
    point[1] = 0.5;
    point[2] = comm_rank;
    point_ids[0] = comm_rank;
    point_dofs[0] = 0.0;
    points[0] = Point(point_ids[0],comm_rank,point);

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup( 
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
    int num_boxes = 5;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = BoxGeometry(box_ids[i],comm_rank,box_ids[i],
		       0.0,0.0,box_ids[i],1.0,1.0,box_ids[i]+1.0);
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // RANGE SETUP
    // Make a range entity set.
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points+1 );
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
    }

    // Add a bad point.
    point_ids[5] = num_points*comm_rank + 1000;
    point[0] = -100.0;
    point[1] = 0.0;
    point[2] = 0.0;
    points[5] = Point(point_ids[5],comm_rank,point);

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup( 
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
    int num_boxes = 5;
    Teuchos::Array<DataTransferKit::SupportId> box_ids( num_boxes );
    Teuchos::ArrayRCP<double> box_dofs( num_boxes );
    Teuchos::Array<Entity> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	box_ids[i] = num_boxes*(comm_size-comm_rank-1) + i;
	box_dofs[i] = 2.0*box_ids[i];
	boxes[i] = BoxGeometry(box_ids[i],comm_rank,box_ids[i],
		       box_ids[i],box_ids[i],box_ids[i],
		       box_ids[i]+1.0,box_ids[i]+1.0,box_ids[i]+1.0);
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( comm, 3, boxes() );
    
    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  boxes(), 1, box_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // RANGE SETUP
    // Make a range entity set.
    int num_points = 5;
    Teuchos::Array<double> point(3);
    Teuchos::Array<DataTransferKit::SupportId> point_ids( num_points+1 );
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
    }

    // Add a bad point.
    DataTransferKit::SupportId id = num_points*comm_rank;
    point_ids[5] = id + 1000;
    point[0] = id + 0.5;
    point[1] = id + 0.5;
    point[2] = id + 1.5;
    points[5] = Point(point_ids[5],comm_rank,point);

    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( comm, 3, points() );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  points(), 1, point_dofs,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_dofs =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // MAPPING
    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    parameters->sublist("Consistent Interpolation");
    Teuchos::ParameterList& search_list = parameters->sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",true);
    Teuchos::RCP<ConsistentInterpolationOperator> map_op = Teuchos::rcp(
	new ConsistentInterpolationOperator(
	    domain_dofs->getMap(),range_dofs->getMap(),*parameters) );

    // Setup the map.
    map_op->setup( 
	domain_manager.functionSpace(), range_manager.functionSpace() );

    // Apply the map.
    map_op->apply( *domain_dofs, *range_dofs );

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
