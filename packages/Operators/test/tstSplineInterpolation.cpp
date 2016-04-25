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
 * \file   tstSplineInterpolation.cpp
 * \author Stuart R. Slattery
 * \brief  SplineInterpolation tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <cstdlib>

#include <DTK_SplineInterpolationOperator.hpp>
#include <DTK_WuBasis.hpp>
#include <DTK_Point.hpp>
#include <DTK_BasicGeometryManager.hpp>
#include <DTK_Entity.hpp>
#include <DTK_EntityCenteredField.hpp>
#include <DTK_FieldMultiVector.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "BelosTypes.hpp"

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationOperator, spline_radius_test )
{
    // Get the test parameters.
    Teuchos::RCP<Teuchos::ParameterList> parameters =
	Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile(
	"spline_interpolation_test_radius.xml", Teuchos::inoutArg(*parameters) );

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    int inverse_rank = comm_size - comm_rank - 1;
    const int space_dim = 3;
    int field_dim = 1;

    // Seed the RNG
    std::srand( 324903231 );

    // Make a set of domain points.
    int num_points = 10;
    int domain_mult = 100;
    int num_domain_points = num_points * domain_mult;
    Teuchos::Array<DataTransferKit::Entity> domain_points( num_domain_points );
    Teuchos::Array<double> coords( space_dim );
    DataTransferKit::EntityId point_id = 0;
    Teuchos::ArrayRCP<double> domain_data( field_dim*num_domain_points );
    for ( int i = 0; i < num_domain_points; ++i )
    {
	point_id = num_domain_points*comm_rank + i;
	coords[0] = (double) std::rand() / (double) RAND_MAX + comm_rank;
	coords[1] = (double) std::rand() / (double) RAND_MAX;
	coords[2] = (double) std::rand() / (double) RAND_MAX;
	domain_points[i] = DataTransferKit::Point( point_id, comm_rank, coords );
	domain_data[i] = coords[0];
    }

    // Make a set of range points.
    Teuchos::Array<DataTransferKit::Entity> range_points( num_points );
    Teuchos::ArrayRCP<double> range_data( field_dim*num_points );
    Teuchos::Array<double> gold_data( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_id = num_points*inverse_rank + i + 1;
	coords[0] = (double) std::rand() / (double) RAND_MAX + inverse_rank;
	coords[1] = (double) std::rand() / (double) RAND_MAX;
	coords[2] = (double) std::rand() / (double) RAND_MAX;
	range_points[i] = DataTransferKit::Point( point_id, comm_rank, coords );
	range_data[i] = 0.0;
	gold_data[i] = coords[0];
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( 
	comm, space_dim, domain_points() );
					   
    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( 
	comm, space_dim, range_points() );

    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  domain_points(), field_dim, domain_data,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_vector =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  range_points(), field_dim, range_data,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_vector =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // Make a spline interpolation operator.
    Teuchos::RCP<DataTransferKit::MapOperator> spline_op =
	Teuchos::rcp( new DataTransferKit::SplineInterpolationOperator<
		      DataTransferKit::WuBasis<2>,space_dim
		      >(domain_vector->getMap(),range_vector->getMap(),*parameters) );

    // Setup the operator.
    spline_op->setup( domain_manager.functionSpace(),
		      range_manager.functionSpace() );

    // Apply the operator.
    spline_op->apply( *domain_vector, *range_vector );

    // Check the apply.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_FLOATING_EQUALITY( gold_data[i], range_data[i], 1.0e-6 );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationOperator, spline_knn_test )
{
    // Get the test parameters.
    Teuchos::RCP<Teuchos::ParameterList> parameters =
	Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile(
	"spline_interpolation_test_knn.xml", Teuchos::inoutArg(*parameters) );

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    int inverse_rank = comm_size - comm_rank - 1;
    const int space_dim = 3;
    int field_dim = 1;

    // Seed the RNG
    std::srand( 324903231 );

    // Make a set of domain points.
    int num_points = 10;
    int domain_mult = 100;
    int num_domain_points = num_points * domain_mult;
    Teuchos::Array<DataTransferKit::Entity> domain_points( num_domain_points );
    Teuchos::Array<double> coords( space_dim );
    DataTransferKit::EntityId point_id = 0;
    Teuchos::ArrayRCP<double> domain_data( field_dim*num_domain_points );
    for ( int i = 0; i < num_domain_points; ++i )
    {
	point_id = num_domain_points*comm_rank + i;
	coords[0] = (double) std::rand() / (double) RAND_MAX + comm_rank;
	coords[1] = (double) std::rand() / (double) RAND_MAX;
	coords[2] = (double) std::rand() / (double) RAND_MAX;
	domain_points[i] = DataTransferKit::Point( point_id, comm_rank, coords );
	domain_data[i] = coords[0];
    }

    // Make a set of range points.
    Teuchos::Array<DataTransferKit::Entity> range_points( num_points );
    Teuchos::ArrayRCP<double> range_data( field_dim*num_points );
    Teuchos::Array<double> gold_data( num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_id = num_points*inverse_rank + i + 1;
	coords[0] = (double) std::rand() / (double) RAND_MAX + inverse_rank;
	coords[1] = (double) std::rand() / (double) RAND_MAX;
	coords[2] = (double) std::rand() / (double) RAND_MAX;
	range_points[i] = DataTransferKit::Point( point_id, comm_rank, coords );
	range_data[i] = 0.0;
	gold_data[i] = coords[0];
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( 
	comm, space_dim, domain_points() );
					   
    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( 
	comm, space_dim, range_points() );

    // Make a DOF vector for the domain.
    Teuchos::RCP<DataTransferKit::Field> domain_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  domain_points(), field_dim, domain_data,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > domain_vector =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  domain_field,
			  domain_manager.functionSpace()->entitySet()) );

    // Make a DOF vector for the range.
    Teuchos::RCP<DataTransferKit::Field> range_field =
	Teuchos::rcp( new DataTransferKit::EntityCenteredField(
			  range_points(), field_dim, range_data,
			  DataTransferKit::EntityCenteredField::BLOCKED) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > range_vector =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  range_field,
			  range_manager.functionSpace()->entitySet()) );

    // Make a spline interpolation operator.
    Teuchos::RCP<DataTransferKit::MapOperator> spline_op =
	Teuchos::rcp( new DataTransferKit::SplineInterpolationOperator<
		      DataTransferKit::WuBasis<2>,space_dim
		      >(domain_vector->getMap(),range_vector->getMap(),*parameters) );

    // Setup the operator.
    spline_op->setup( domain_manager.functionSpace(),
		      range_manager.functionSpace() );

    // Apply the operator.
    spline_op->apply( *domain_vector, *range_vector );

    // Check the apply.
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_FLOATING_EQUALITY( gold_data[i], range_data[i], 1.0e-6 );
    }
}

//---------------------------------------------------------------------------//
// end tstSplineInterpolation.cpp
//---------------------------------------------------------------------------//
