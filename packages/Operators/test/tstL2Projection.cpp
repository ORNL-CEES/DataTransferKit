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
 * \file tstL2ProjectionOperator.cpp
 * \author Stuart R. Slattery
 * \brief L2ProjectionOperator unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "reference_implementation/DTK_ReferenceHexMesh.hpp"

#include <DTK_FieldMultiVector.hpp>

#include <DTK_L2ProjectionOperator.hpp>

#include <DTK_BasicEntityPredicates.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_Import.hpp>

//---------------------------------------------------------------------------//
// TEST EPSILON
//---------------------------------------------------------------------------//

// Floating point epsilon for checking the transferred field.
const double field_epsilon = 1.0e-7;

// Floating point epsilon for checking the integrated field.
const double integral_epsilon = 1.0e-10;

//---------------------------------------------------------------------------//
// TEST FUNCTION
//---------------------------------------------------------------------------//
double testFunction( const Teuchos::ArrayView<double> &coords )
{
    return 9.3 * coords[0] + 2.2 * coords[1] + 1.33 * coords[2] + 1.0;
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS.
//---------------------------------------------------------------------------//
// integrate a field over the mesh.
double integrateField( DataTransferKit::UnitTest::ReferenceHexMesh &mesh,
                       DataTransferKit::Field &field )
{
    // Get the entity set.
    auto set = mesh.functionSpace()->entitySet();

    // Get the communicator.
    auto comm = set->communicator();

    // Import the field to a ghosted decomposition so we have access to node
    // DOFs that are not locally owned.
    Teuchos::RCP<DataTransferKit::FieldMultiVector> vector =
        Teuchos::rcp( new DataTransferKit::FieldMultiVector(
            comm, Teuchos::rcpFromRef( field ) ) );
    vector->pullDataFromApplication();

    auto ghosted_field = mesh.ghostedNodalField( field.dimension() );
    Teuchos::RCP<DataTransferKit::FieldMultiVector> ghosted_vector =
        Teuchos::rcp(
            new DataTransferKit::FieldMultiVector( comm, ghosted_field ) );
    Tpetra::Import<DataTransferKit::FieldMultiVector::LO,
                   DataTransferKit::FieldMultiVector::GO,
                   DataTransferKit::FieldMultiVector::Node>
        importer( vector->getMap(), ghosted_vector->getMap() );
    ghosted_vector->doImport( *vector, importer, Tpetra::REPLACE );
    ghosted_vector->pushDataToApplication();

    // Get the cells.
    int space_dim = 3;
    DataTransferKit::LocalEntityPredicate local_pred( comm->getRank() );
    auto cells = set->entityIterator( space_dim, local_pred.getFunction() );

    // Get the integration weights and points. All entities are a linear hex
    // so we can cache this here using the first cell.
    auto ir = mesh.functionSpace()->integrationRule();
    Teuchos::Array<Teuchos::Array<double>> points;
    Teuchos::Array<double> weights;
    int integration_order = 3;
    ir->getIntegrationRule( *cells, integration_order, points, weights );
    int num_points = weights.size();

    // Get the shape function.
    auto shape = mesh.functionSpace()->shapeFunction();
    Teuchos::Array<DataTransferKit::SupportId> dofs;
    Teuchos::Array<double> shape_evals;

    // Get the local map.
    auto local_map = mesh.functionSpace()->localMap();

    // Calculate the local integral.
    int num_nodes = 8;
    double local_integral = 0.0;
    auto cells_begin = cells.begin();
    auto cells_end = cells.end();
    for ( cells = cells_begin; cells != cells_end; ++cells )
    {
        // Get the cell dofs.
        shape->entitySupportIds( *cells, dofs );

        // Add the cell contribution to the integral.
        double cell_integral = 0.0;
        for ( int p = 0; p < num_points; ++p )
        {
            // Evaluate the basis at the integration point.
            shape->evaluateValue( *cells, points[p](), shape_evals );

            // Compute the field value at the integration point.
            double point_value = 0.0;
            for ( int n = 0; n < num_nodes; ++n )
            {
                point_value +=
                    shape_evals[n] * ghosted_field->readFieldData( dofs[n], 0 );
            }

            // Add to the cell integral.
            cell_integral += weights[p] * point_value;
        }

        // Add the cell integral scaled by the cell volume.
        local_integral += cell_integral * local_map->measure( *cells );
    }

    // Compute the global integral.
    double global_integral = 0.0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_integral,
                        Teuchos::ptrFromRef( global_integral ) );

    // Return the integral. We divide by 8 because of how intrepid scales the
    // integration weights.
    return global_integral / 8.0;
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( L2ProjectionOperator, l2_projection )
{
    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    DataTransferKit::LocalEntityPredicate local_pred( comm->getRank() );

    // Set the global problem bounds.
    double x_min = 0.0;
    double y_min = 0.0;
    double z_min = 0.0;
    double x_max = 3.1;
    double y_max = 5.2;
    double z_max = 8.3;

    // Create a target mesh and field.
    int num_sx = 8;
    int num_sy = 8;
    int num_sz = 8;

    DataTransferKit::UnitTest::ReferenceHexMesh source_mesh(
        comm, x_min, x_max, num_sx, y_min, y_max, num_sy, z_min, z_max,
        num_sz );
    auto source_field = source_mesh.nodalField( 1 );
    Teuchos::RCP<DataTransferKit::FieldMultiVector> source_vector =
        Teuchos::rcp(
            new DataTransferKit::FieldMultiVector( comm, source_field ) );

    // Put some data on the source field.
    auto source_local_map = source_mesh.functionSpace()->localMap();
    auto source_nodes =
        source_mesh.functionSpace()->entitySet()->entityIterator(
            0, local_pred.getFunction() );
    auto source_nodes_begin = source_nodes.begin();
    auto source_nodes_end = source_nodes.end();
    Teuchos::Array<double> source_coords( 3 );
    for ( source_nodes = source_nodes.begin();
          source_nodes != source_nodes.end(); ++source_nodes )
    {
        int i, j, k;
        source_mesh.id( source_nodes->id(), i, j, k );

        source_local_map->centroid( *source_nodes, source_coords() );

        source_field->writeFieldData( source_nodes->id(), 0,
                                      testFunction( source_coords() ) );
    }

    // Create a source mesh and field.
    int num_tx = 9;
    int num_ty = 7;
    int num_tz = 7;
    DataTransferKit::UnitTest::ReferenceHexMesh target_mesh(
        comm, x_min, x_max, num_tx, y_min, y_max, num_ty, z_min, z_max,
        num_tz );
    auto target_field = target_mesh.nodalField( 1 );
    Teuchos::RCP<DataTransferKit::FieldMultiVector> target_vector =
        Teuchos::rcp(
            new DataTransferKit::FieldMultiVector( comm, target_field ) );

    // Create a map.
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::parameterList();
    Teuchos::ParameterList &l2_list = parameters->sublist( "L2 Projection" );
    l2_list.set( "Integration Order", 3 );
    Teuchos::ParameterList &search_list = parameters->sublist( "Search" );
    search_list.set( "Point Inclusion Tolerance", 1.0e-6 );

    Teuchos::RCP<DataTransferKit::L2ProjectionOperator> map_op =
        Teuchos::rcp( new DataTransferKit::L2ProjectionOperator(
            source_vector->getMap(), target_vector->getMap(), *parameters ) );

    // Setup the map.
    map_op->setup( source_mesh.functionSpace(), target_mesh.functionSpace() );

    // Apply the map.
    map_op->apply( *source_vector, *target_vector );

    // Check the results of the mapping.
    auto target_nodes =
        target_mesh.functionSpace()->entitySet()->entityIterator(
            0, local_pred.getFunction() );
    auto target_nodes_begin = target_nodes.begin();
    auto target_nodes_end = target_nodes.end();
    auto target_local_map = target_mesh.functionSpace()->localMap();
    Teuchos::Array<double> target_coords( 3 );
    for ( target_nodes = target_nodes.begin();
          target_nodes != target_nodes.end(); ++target_nodes )
    {
        int i, j, k;
        target_mesh.id( target_nodes->id(), i, j, k );

        target_local_map->centroid( *target_nodes, target_coords() );
        double gold_data = testFunction( target_coords() );
        double target_data =
            target_field->readFieldData( target_nodes->id(), 0 );
        TEST_FLOATING_EQUALITY( target_data, gold_data, field_epsilon );
    }

    // Check accurate preservation of the global integral. We are using the
    // source mesh quadrature so this should be very accurate.
    double source_integral = integrateField( source_mesh, *source_field );
    double target_integral = integrateField( target_mesh, *target_field );
    TEST_FLOATING_EQUALITY( source_integral, target_integral,
                            integral_epsilon );
}

//---------------------------------------------------------------------------//
// end tstL2ProjectionOperator.cpp
//---------------------------------------------------------------------------//
