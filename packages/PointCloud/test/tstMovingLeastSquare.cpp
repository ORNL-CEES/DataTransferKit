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
 * \file   tstMovingLeastSquare.cpp
 * \author Stuart R. Slattery
 * \brief  MovingLeastSquare tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <limits>

#include <DTK_MovingLeastSquareReconstructionOperator.hpp>
#include <DTK_WuBasis.hpp>
#include <DTK_Point.hpp>
#include <DTK_BasicGeometryManager.hpp>
#include <DTK_Entity.hpp>
#include <DTK_EntityCenteredDOFVector.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MovingLeastSquareReconstructionOperator, mls_test )
{
    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    int inverse_rank = comm_size - comm_rank - 1;
    const int space_dim = 3;
    int field_dim = 1;

    // Make a set of domain points.
    int num_points = 10;
    Teuchos::Array<DataTransferKit::Entity> domain_points( num_points );
    Teuchos::Array<double> coords( space_dim );
    DataTransferKit::EntityId point_id = 0;
    Teuchos::ArrayRCP<double> domain_data( field_dim*num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_id = num_points*comm_rank + i;
	coords[0] = point_id;
	coords[1] = point_id;
	coords[2] = point_id;
	domain_points[i] = DataTransferKit::Point( point_id, comm_rank, coords );
	domain_data[i] = point_id;
    }

    // Make a set of range points.
    Teuchos::Array<DataTransferKit::Entity> range_points( num_points );
    Teuchos::ArrayRCP<double> range_data( field_dim*num_points );
    for ( int i = 0; i < num_points; ++i )
    {
	point_id = num_points*inverse_rank + i;
	coords[0] = point_id;
	coords[1] = point_id;
	coords[2] = point_id;
	range_points[i] = DataTransferKit::Point( point_id, comm_rank, coords );
	range_data[i] = 0.0;
    }

    // Make a manager for the domain geometry.
    DataTransferKit::BasicGeometryManager domain_manager( 
	comm, space_dim, DataTransferKit::ENTITY_TYPE_NODE, domain_points() );
					   
    // Make a manager for the range geometry.
    DataTransferKit::BasicGeometryManager range_manager( 
	comm, space_dim, DataTransferKit::ENTITY_TYPE_NODE, range_points() );

    // Make a DOF vector for the domain.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > domain_vector =
	DataTransferKit::EntityCenteredDOFVector::createTpetraMultiVectorFromEntitiesAndView(
	    comm, domain_points(), field_dim, domain_data );

    // Make a DOF vector for the range.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > range_vector =
	DataTransferKit::EntityCenteredDOFVector::createTpetraMultiVectorFromEntitiesAndView(
	    comm, range_points(), field_dim, range_data );

    // Make a moving least square reconstruction operator.
    double radius = 0.25;
    Teuchos::RCP<DataTransferKit::MapOperator<double> > mls_op =
	Teuchos::rcp( 
	    new DataTransferKit::MovingLeastSquareReconstructionOperator<double,DataTransferKit::WuBasis<2>,space_dim>(radius) );

    // Setup the operator.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    mls_op->setup( domain_vector->getMap(),
		   domain_manager.functionSpace(),
		   range_vector->getMap(),
		   range_manager.functionSpace(),
		   parameters );

    // Apply the operator.
    mls_op->apply( *domain_vector, *range_vector );

    // Check the apply.
    for ( int i = 0; i < num_points; ++i )
    {
	double test_val = num_points*inverse_rank + i;
	TEST_FLOATING_EQUALITY( range_data[i], test_val,
				10.0*std::numeric_limits<double>::epsilon() );
    }
}

//---------------------------------------------------------------------------//
// end tstMovingLeastSquare.cpp
//---------------------------------------------------------------------------//
