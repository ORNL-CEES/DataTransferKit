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
 * \file tstIntegral.cpp
 * \author Stuart R. Slattery
 * \brief Integrator unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_Import.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>

#include "reference_implementation/DTK_ReferenceHexMesh.hpp"

#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_DBC.hpp"
#include "DTK_IntegrationPoint.hpp"
#include "DTK_Jacobian.hpp"

#define NX 8 // number of cells in x direction
#define NY 8 // number of cells in y direction
#define NZ 8 // number of cells in z direction

//---------------------------------------------------------------------------//
// TEST FUNCTION
//---------------------------------------------------------------------------//
int FUNCTION = 0;
int NUM_FUNCTIONS = 4;
double testFunction( const Teuchos::ArrayView<double> &coords )
{
    DTK_REQUIRE( coords.size() == 3 );
    DTK_REQUIRE( FUNCTION >= 0 && FUNCTION < NUM_FUNCTIONS );

    double x = coords[0], y = coords[1], z = coords[2];

    if ( FUNCTION == 0 )
        return 1.0;
    else if ( FUNCTION == 1 )
        return 9.3 * x + 2.2 * y + 3.4 * z + 1.0;
    else if ( FUNCTION == 2 )
        return x * y;
    else if ( FUNCTION == 3 )
    {
        // align with non-perturbed grid to validate
        double a = 3.1 / 2;
        return std::abs( x - a );
    }

    // should not reach this
    DTK_CHECK( false );
}
double integralExact( double x_min, double x_max, double y_min, double y_max,
                      double z_min, double z_max )
{
    DTK_REQUIRE( FUNCTION >= 0 && FUNCTION < NUM_FUNCTIONS );

    if ( FUNCTION == 0 )
        return ( x_max - x_min ) * ( y_max - y_min ) * ( z_max - z_min );
    else if ( FUNCTION == 1 )
        return 9.3 / 2 * ( x_max * x_max - x_min * x_min ) * ( y_max - y_min ) *
                   ( z_max - z_min ) +
               2.2 / 2 * ( y_max * y_max - y_min * y_min ) * ( x_max - x_min ) *
                   ( z_max - z_min ) +
               3.4 / 2 * ( z_max * z_max - z_min * z_min ) * ( x_max - x_min ) *
                   ( y_max - y_min ) +
               ( x_max - x_min ) * ( y_max - y_min ) * ( z_max - z_min );
    else if ( FUNCTION == 2 )
        return 1.0 / 4 * ( x_max * x_max - x_min * x_min ) *
               ( y_max * y_max - y_min * y_min ) * ( z_max - z_min );
    else if ( FUNCTION == 3 )
    {
        DTK_REQUIRE( x_min < 2 && x_max > 2 );
        double a = 3.1 / 2;
        return 0.5 * ( ( x_min - a ) * ( x_min - a ) +
                       ( x_max - a ) * ( x_max - a ) ) *
               ( y_max - y_min ) * ( z_max - z_min );
    }

    // should not reach this
    DTK_CHECK( false );
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS.
//---------------------------------------------------------------------------//
// Given an element id compute the ids of the associated nodes.
void compute_elem_node_ids( int elem_id, int node_ids[] )
{
    int I = elem_id % NX, i = I;
    int J = ( elem_id / NX ) % NY, j = J;
    int K = elem_id / ( NX * NY ), k = K;
    int nx = NX + 1, ny = NY + 1;

    int base = k * nx * ny + j * nx + i;
    node_ids[0] = base;
    node_ids[1] = base + 1;
    node_ids[2] = base + 1 + nx;
    node_ids[3] = base + nx;
    node_ids[4] = base + nx * ny;
    node_ids[5] = base + nx * ny + 1;
    node_ids[6] = base + nx * ny + 1 + nx;
    node_ids[7] = base + nx * ny + nx;
}

enum DTKVersion
{
    // The incorrect way to integrate a function over the mesh (the original L2
    // approach)
    DTK_OLD,
    // The incorrect way to integrate a function over the mesh but using a
    // correct measure.
    DTK_TRANSITION,
    // The correct way to integrate a function over the mesh, using weights
    // multiplied by Jacobians
    DTK_NEW
};

double integrateDTK( const Teuchos::RCP<DataTransferKit::FunctionSpace> &space,
                     int integration_order, DTKVersion version )
{
    // Get the cells.
    int space_dim = 3;

    // Get the data.
    auto set = space->entitySet();
    auto shape_function = space->shapeFunction();
    auto local_map = space->localMap();
    auto integration_rule = space->integrationRule();
    auto cells = space->entitySet()->entityIterator( space_dim );

    DataTransferKit::Jacobian jacobian( space );

    // Calculate the local integral.
    double integral = 0.0;

    auto cells_begin = cells.begin();
    auto cells_end = cells.end();
    for ( auto it = cells_begin; it != cells_end; ++it )
    {
        // Get the support ids of the entity.
        Teuchos::Array<DataTransferKit::SupportId> support_ids;
        shape_function->entitySupportIds( *it, support_ids );

        // Get the integration rule.
        Teuchos::Array<Teuchos::Array<double>> int_points;
        Teuchos::Array<double> int_weights;
        integration_rule->getIntegrationRule( *it, integration_order,
                                              int_points, int_weights );
        int num_ip = int_weights.size();

        // Compute the measure of the entity.
        double entity_measure = 0.0;
        if ( version == DTK_OLD )
        {
            entity_measure = local_map->measure( *it );
        }
        else if ( version == DTK_TRANSITION || version == DTK_NEW )
        {
            for ( int p = 0; p < num_ip; ++p )
            {
                // Compute derminant.
                double det =
                    jacobian.jacobian_determinant( *it, int_points[p]() );

                // Update the measure.
                entity_measure += int_weights[p] * det;
            }
        }

        // Create new integration points.
        for ( int p = 0; p < num_ip; ++p )
        {
            DataTransferKit::IntegrationPoint ip;
            ip.d_physical_coordinates.resize( space_dim );
            ip.d_owner_measure = entity_measure;
            ip.d_owner_support_ids = support_ids;
            ip.d_integration_weight = int_weights[p];

            // Evaluate the shape function.
            Teuchos::Array<double> shape_evals;
            shape_function->evaluateValue( *it, int_points[p](), shape_evals );
            ip.d_owner_shape_evals = shape_evals;

            // Map the integration point to the physical frame of the range
            // entity.
            local_map->mapToPhysicalFrame( *it, int_points[p](),
                                           ip.d_physical_coordinates() );

            // Update the integral.
            if ( version == DTK_OLD || version == DTK_TRANSITION )
            {
                integral += entity_measure * int_weights[p] *
                            testFunction( ip.d_physical_coordinates() );
            }
            else if ( version == DTK_NEW )
            {
                // Compute derminant.
                double det =
                    jacobian.jacobian_determinant( *it, int_points[p]() );
                integral += int_weights[p] * det *
                            testFunction( ip.d_physical_coordinates() );
            }
        }
    }

    if ( version == DTK_OLD || version == DTK_TRANSITION )
    {
        integral /= 8.0;
    }

    // Return the integral.
    return integral;
}

double
integrateDTKField( const Teuchos::RCP<DataTransferKit::FunctionSpace> &space,
                   const Teuchos::RCP<DataTransferKit::Field> &field,
                   int integration_order, DTKVersion version )
{
    // Get the cells.
    int space_dim = 3;

    // Get the data.
    auto set = space->entitySet();
    auto shape_function = space->shapeFunction();
    auto local_map = space->localMap();
    auto integration_rule = space->integrationRule();
    auto cells = space->entitySet()->entityIterator( space_dim );

    DataTransferKit::Jacobian jacobian( space );

    // Calculate the local integral.
    double integral = 0.0;

    auto cells_begin = cells.begin();
    auto cells_end = cells.end();
    for ( auto it = cells_begin; it != cells_end; ++it )
    {
        // Get the support ids of the entity.
        Teuchos::Array<DataTransferKit::SupportId> support_ids;
        shape_function->entitySupportIds( *it, support_ids );
        int cardinality = support_ids.size();

        // Get the integration rule.
        Teuchos::Array<Teuchos::Array<double>> int_points;
        Teuchos::Array<double> int_weights;
        integration_rule->getIntegrationRule( *it, integration_order,
                                              int_points, int_weights );
        int num_ip = int_weights.size();

        // Compute the measure of the entity.
        double entity_measure = 0.0;
        if ( version == DTK_OLD )
        {
            entity_measure = local_map->measure( *it );
        }
        else if ( version == DTK_TRANSITION || version == DTK_NEW )
        {
            for ( int p = 0; p < num_ip; ++p )
            {
                // Compute derminant.
                double det =
                    jacobian.jacobian_determinant( *it, int_points[p]() );

                // Update the measure.
                entity_measure += int_weights[p] * det;
            }
        }

        // Create new integration points.
        for ( int p = 0; p < num_ip; ++p )
        {
            // Evaluate the shape function.
            Teuchos::Array<double> shape_evals;
            shape_function->evaluateValue( *it, int_points[p](), shape_evals );

            // Compute the field value at the integration point.
            double point_value = 0.0;
            for ( int n = 0; n < cardinality; ++n )
            {
                point_value +=
                    shape_evals[n] * field->readFieldData( support_ids[n], 0 );
            }

            // Update the integral.
            if ( version == DTK_OLD || version == DTK_TRANSITION )
            {
                integral += entity_measure * int_weights[p] * point_value;
            }
            else if ( version == DTK_NEW )
            {
                // Compute derminant.
                double det =
                    jacobian.jacobian_determinant( *it, int_points[p]() );
                integral += int_weights[p] * det * point_value;
            }
        }
    }

    if ( version == DTK_OLD || version == DTK_TRANSITION )
    {
        integral /= 8.0;
    }

    // Return the integral.
    return integral;
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Integrator, integration )
{
    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();

    DTK_REQUIRE( comm->getSize() == 1 );

    // Set the global problem bounds.
    double x_min = 0.0;
    double y_min = 0.0;
    double z_min = 0.0;
    double x_max = 3.1;
    double y_max = 5.2;
    double z_max = 8.3;

    // Create a target mesh.
    int num_sx = NX;
    int num_sy = NY;
    int num_sz = NZ;

    // Prin the header.
    for ( FUNCTION = 0; FUNCTION < NUM_FUNCTIONS; FUNCTION++ )
    {
        std::cout << std::endl;
        switch ( FUNCTION )
        {
        case 0:
            std::cout << "*f(x) = 1*" << std::endl;
            break;
        case 1:
            std::cout << "*f(x) = 9.3 x + 2.2 y + 3.4 z + 1*" << std::endl;
            break;
        case 2:
            std::cout << "*f(x) = x y*" << std::endl;
            break;
        case 3:
            std::cout << "*f(x) = |x - 1.55|*" << std::endl;
        }

        for ( int integration_order = 1; integration_order <= 2;
              integration_order++ )
        {
            std::cout << "\nintegration order = " << integration_order
                      << std::endl;

            std::cout << std::endl;
            std::cout << "pert |  DTK old  |  DTK trs  |  DTK new"
                      << "  | fDTK old  | fDTK trs  | fDTK new" << std::endl;
            std::cout << "---- | --------- | --------- | ---------"
                      << " | --------- | --------- | ---------" << std::endl;
            for ( double perturb = 0.0; perturb < 0.45; perturb += 0.1 )
            {
                // Create the mesh.
                DataTransferKit::UnitTest::ReferenceHexMesh mesh(
                    comm, x_min, x_max, num_sx, y_min, y_max, num_sy, z_min,
                    z_max, num_sz, perturb );
                auto space = mesh.functionSpace();

                // Create the field.
                auto field = mesh.nodalField( 1 );

                // Put test data on the field.
                auto local_map = mesh.functionSpace()->localMap();
                auto nodes = space->entitySet()->entityIterator( 0 );

                Teuchos::Array<double> coords( 3 );
                for ( auto node = nodes.begin(); node != nodes.end(); ++node )
                {
                    int i, j, k;
                    mesh.id( node->id(), i, j, k );

                    local_map->centroid( *node, coords() );

                    field->writeFieldData( node->id(), 0,
                                           testFunction( coords() ) );
                }

                // Compare two integrals
                double integral_exact =
                    integralExact( x_min, x_max, y_min, y_max, z_min, z_max );
                double integral_old =
                    integrateDTK( space, integration_order, DTK_OLD );
                double integral_trs =
                    integrateDTK( space, integration_order, DTK_TRANSITION );
                double integral_new =
                    integrateDTK( space, integration_order, DTK_NEW );
                double integral_old_f = integrateDTKField(
                    space, field, integration_order, DTK_OLD );
                double integral_trs_f = integrateDTKField(
                    space, field, integration_order, DTK_TRANSITION );
                double integral_new_f = integrateDTKField(
                    space, field, integration_order, DTK_NEW );

                std::cout << std::fixed << std::setprecision( 2 ) << perturb
                          << " | " << std::scientific << std::setprecision( 3 )
                          << std::abs( integral_old - integral_exact ) << " | "
                          << std::abs( integral_trs - integral_exact ) << " | "
                          << std::abs( integral_new - integral_exact ) << " | "
                          << std::abs( integral_old_f - integral_exact )
                          << " | "
                          << std::abs( integral_trs_f - integral_exact )
                          << " | "
                          << std::abs( integral_new_f - integral_exact )
                          << std::endl;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// end tstIntegral.cpp
//---------------------------------------------------------------------------//
