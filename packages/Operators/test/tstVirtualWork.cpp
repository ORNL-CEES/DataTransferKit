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
 * \file   tstVirtualWork.cpp
 * \author Stuart R. Slattery
 * \brief  Virtual work property tests for point cloud operators.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <DTK_BasicGeometryManager.hpp>
#include <DTK_Entity.hpp>
#include <DTK_EntityCenteredField.hpp>
#include <DTK_FieldMultiVector.hpp>
#include <DTK_MapOperatorFactory.hpp>
#include <DTK_Point.hpp>

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

//---------------------------------------------------------------------------//
// Test epsilon.
//---------------------------------------------------------------------------//

const double epsilon = 1.0e-14;

//---------------------------------------------------------------------------//
// Test dirver.
//---------------------------------------------------------------------------//
// This test checks the concept of virtual work conservation in fluid
// structure interaction problems. The point cloud operators are designed to
// preserve virtual work. Virtual work is effectively the sum of all the force
// and displacement products at each node on the coupled surface (i.e. virtual
// work = displacement * force). This quantity is dimension agnostic so we
// check it here in three dimensions.
void setupAndRunTest( const std::string &input_file,
                      double &structure_virtual_work,
                      double &fluid_virtual_work )
{
    // Get the test parameters.
    Teuchos::RCP<Teuchos::ParameterList> parameters =
        Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( input_file,
                                          Teuchos::inoutArg( *parameters ) );

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    int inverse_rank = comm_size - comm_rank - 1;
    const int space_dim = 3;

    // Seed the RNG
    std::srand( 324903231 );

    // Make a set of fluid points. These span 0-1 in y, and z and span
    // comm_rank-comm_rank+1 in x. The value of the field we are transferring
    // is the x + y + z coordinate of the points.
    int num_points = 1000;
    Teuchos::Array<DataTransferKit::Entity> fluid_points( num_points );
    Teuchos::Array<double> coords( space_dim );
    DataTransferKit::EntityId point_id = 0;
    Teuchos::ArrayRCP<double> fluid_displacements( space_dim * num_points );
    Teuchos::ArrayRCP<double> fluid_forces( space_dim * num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        point_id = num_points * comm_rank + i;
        coords[0] = (double)std::rand() / (double)RAND_MAX + comm_rank;
        coords[1] = (double)std::rand() / (double)RAND_MAX;
        coords[2] = (double)std::rand() / (double)RAND_MAX;
        fluid_points[i] = DataTransferKit::Point( point_id, comm_rank, coords );

        for ( int d = 0; d < space_dim; ++d )
        {
            fluid_displacements[num_points * d + i] = 0.0;
            fluid_forces[num_points * d + i] = coords[d] + 1.5 * d;
        }
    }

    // Make a set of structure points. These span 0-1 in y and z and span
    // comm_rank-inverse_rank+1 in x.
    Teuchos::Array<DataTransferKit::Entity> structure_points( num_points );
    Teuchos::ArrayRCP<double> structure_displacements( space_dim * num_points );
    Teuchos::ArrayRCP<double> structure_forces( space_dim * num_points );
    for ( int i = 0; i < num_points; ++i )
    {
        point_id = num_points * inverse_rank + i;
        coords[0] = (double)std::rand() / (double)RAND_MAX + inverse_rank;
        coords[1] = (double)std::rand() / (double)RAND_MAX;
        coords[2] = (double)std::rand() / (double)RAND_MAX;
        structure_points[i] =
            DataTransferKit::Point( point_id, comm_rank, coords );

        for ( int d = 0; d < space_dim; ++d )
        {
            structure_displacements[num_points * d + i] =
                coords[space_dim - d - 1] + 2.0 * d;
            structure_forces[num_points * d + i] = 0.0;
        }
    }

    // Make a manager for the fluid geometry.
    DataTransferKit::BasicGeometryManager fluid_manager( comm, space_dim,
                                                         fluid_points() );

    // Make a manager for the structure geometry.
    DataTransferKit::BasicGeometryManager structure_manager(
        comm, space_dim, structure_points() );

    // Make DOF vectors for the fluid fields.
    Teuchos::RCP<DataTransferKit::Field> fluid_displacements_field =
        Teuchos::rcp( new DataTransferKit::EntityCenteredField(
            fluid_points(), space_dim, fluid_displacements,
            DataTransferKit::EntityCenteredField::BLOCKED ) );
    auto fluid_displacements_vector =
        Teuchos::rcp( new DataTransferKit::FieldMultiVector(
            fluid_displacements_field,
            fluid_manager.functionSpace()->entitySet() ) );

    Teuchos::RCP<DataTransferKit::Field> fluid_forces_field =
        Teuchos::rcp( new DataTransferKit::EntityCenteredField(
            fluid_points(), space_dim, fluid_forces,
            DataTransferKit::EntityCenteredField::BLOCKED ) );
    auto fluid_forces_vector =
        Teuchos::rcp( new DataTransferKit::FieldMultiVector(
            fluid_forces_field, fluid_manager.functionSpace()->entitySet() ) );

    // Make DOF vectors for the structure fields.
    Teuchos::RCP<DataTransferKit::Field> structure_displacements_field =
        Teuchos::rcp( new DataTransferKit::EntityCenteredField(
            structure_points(), space_dim, structure_displacements,
            DataTransferKit::EntityCenteredField::BLOCKED ) );
    auto structure_displacements_vector =
        Teuchos::rcp( new DataTransferKit::FieldMultiVector(
            structure_displacements_field,
            structure_manager.functionSpace()->entitySet() ) );

    Teuchos::RCP<DataTransferKit::Field> structure_forces_field =
        Teuchos::rcp( new DataTransferKit::EntityCenteredField(
            structure_points(), space_dim, structure_forces,
            DataTransferKit::EntityCenteredField::BLOCKED ) );
    auto structure_forces_vector =
        Teuchos::rcp( new DataTransferKit::FieldMultiVector(
            structure_forces_field,
            structure_manager.functionSpace()->entitySet() ) );

    // Create the point cloud operator to map from the structure to the fluid.
    DataTransferKit::MapOperatorFactory factory;
    Teuchos::RCP<DataTransferKit::MapOperator> cloud_op =
        factory.create( structure_displacements_vector->getMap(),
                        fluid_displacements_vector->getMap(), *parameters );

    // Setup the operator.
    cloud_op->setup( structure_manager.functionSpace(),
                     fluid_manager.functionSpace() );

    // Map the displacements.
    cloud_op->apply( *structure_displacements_vector,
                     *fluid_displacements_vector, Teuchos::NO_TRANS );

    // Map the forces.
    cloud_op->apply( *fluid_forces_vector, *structure_forces_vector,
                     Teuchos::TRANS );

    // Calculate structure virtual work.
    double local_structure_virtual_work = 0.0;
    for ( int i = 0; i < num_points; ++i )
    {
        for ( int d = 0; d < space_dim; ++d )
        {
            local_structure_virtual_work +=
                structure_displacements[num_points * d + i] *
                structure_forces[num_points * d + i];
        }
    }
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM,
                        local_structure_virtual_work,
                        Teuchos::ptrFromRef( structure_virtual_work ) );

    // Calculate fluid virtual work.
    double local_fluid_virtual_work = 0.0;
    for ( int i = 0; i < num_points; ++i )
    {
        for ( int d = 0; d < space_dim; ++d )
        {
            local_fluid_virtual_work +=
                fluid_displacements[num_points * d + i] *
                fluid_forces[num_points * d + i];
        }
    }
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_fluid_virtual_work,
                        Teuchos::ptrFromRef( fluid_virtual_work ) );
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VirtualWork, mls_radius_test )
{
    // Run the test.
    double structure_virtual_work = 0.0;
    double fluid_virtual_work = 0.0;
    setupAndRunTest( "mls_test_radius.xml", structure_virtual_work,
                     fluid_virtual_work );

    // Check the results.
    TEST_FLOATING_EQUALITY( fluid_virtual_work, structure_virtual_work,
                            epsilon );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( VirtualWork, mls_knn_test )
{
    // Run the test.
    double structure_virtual_work = 0.0;
    double fluid_virtual_work = 0.0;
    setupAndRunTest( "mls_test_knn.xml", structure_virtual_work,
                     fluid_virtual_work );

    // Check the results.
    TEST_FLOATING_EQUALITY( fluid_virtual_work, structure_virtual_work,
                            epsilon );
}

//---------------------------------------------------------------------------//
// end tstVirtualWork.cpp
//---------------------------------------------------------------------------//
