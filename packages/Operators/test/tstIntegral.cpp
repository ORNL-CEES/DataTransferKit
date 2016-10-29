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

#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_Import.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>

#include "reference_implementation/DTK_ReferenceHexMesh.hpp"

#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_DBC.hpp"
#include "DTK_FieldMultiVector.hpp"
#include "DTK_IntrepidCell.hpp"
#include "DTK_IntegrationPoint.hpp"

#define NX 8 // number of cells in x direction
#define NY 8 // number of cells in y direction
#define NZ 8 // number of cells in z direction

//---------------------------------------------------------------------------//
// TEST FUNCTION
//---------------------------------------------------------------------------//
double testFunction( const Teuchos::ArrayView<double> &coords )
{
    DTK_REQUIRE(coords.size() == 3);
    return 1.0;
    // return 9.3 * coords[0] + 2.2 * coords[1] + 1.33 * coords[2] + 1.0;
    // return 9.3 * coords[0]*coords[2] + 2.2 * coords[1]*coords[0] + 1.33 * coords[2] + 1.0;
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS.
//---------------------------------------------------------------------------//
// Given an element id compute the ids of the associated nodes.
void compute_elem_node_ids( int elem_id, int node_ids[] )
{
    int I = elem_id % NX, i = I;
    int J = (elem_id / NX) % NY, j = J;
    int K = elem_id / (NX*NY), k = K;
    int nx = NX+1, ny = NY+1;

    int base = k*nx*ny + j*nx + i;
    node_ids[0] = base;
    node_ids[1] = base + 1;
    node_ids[2] = base + 1 + nx;
    node_ids[3] = base + nx;
    node_ids[4] = base + nx*ny;
    node_ids[5] = base + nx*ny + 1;
    node_ids[6] = base + nx*ny + 1 + nx;
    node_ids[7] = base + nx*ny + nx;
}

// The incorrect way to integrate a field over the mesh. The code is a modified
// version taken from the L2 projection mass matrix assembly.
double integrateFieldDTKOld( DataTransferKit::UnitTest::ReferenceHexMesh &mesh,
                             DataTransferKit::Field &field, int integration_order )
{
    // Get the cells.
    int space_dim = 3;

    // Get the data.
    auto space            = mesh.functionSpace();
    auto shape_function   = space->shapeFunction();
    auto local_map        = space->localMap();
    auto integration_rule = space->integrationRule();
    auto cells            = space->entitySet()->entityIterator( space_dim );

    // Calculate the local integral.
    double integral  = 0.0;

    auto cells_begin = cells.begin();
    auto cells_end   = cells.end();
    for ( auto it = cells_begin; it != cells_end; ++it )
    {
        // Get the support ids of the entity.
        Teuchos::Array<DataTransferKit::SupportId> support_ids;
        shape_function->entitySupportIds( *it, support_ids );
        int cardinality = support_ids.size();

        // Get the measure of the entity.
        double entity_measure = local_map->measure( *it );

        // Get the integration rule.
        Teuchos::Array<Teuchos::Array<double>> int_points;
        Teuchos::Array<double>                 int_weights;
        integration_rule->getIntegrationRule( *it, integration_order,
                                              int_points, int_weights );
        int num_ip = int_weights.size();

        // Create new integration points.
        for ( int p = 0; p < num_ip; ++p )
        {
            DataTransferKit::IntegrationPoint ip;
            ip.d_physical_coordinates.resize( space_dim );
            ip.d_owner_measure      = entity_measure;
            ip.d_owner_support_ids  = support_ids;
            ip.d_integration_weight = int_weights[p];

            // Evaluate the shape function.
            Teuchos::Array<double> shape_evals;
            shape_function->evaluateValue( *it, int_points[p](),
                                           shape_evals );
            ip.d_owner_shape_evals = shape_evals;

            // Map the integration point to the physical frame of the range
            // entity.
            local_map->mapToPhysicalFrame(
                *it, int_points[p](), ip.d_physical_coordinates() );

            // Update the integral.
            for ( int ni = 0; ni < cardinality; ++ni )
                integral += entity_measure * int_weights[p] * testFunction(ip.d_physical_coordinates());
        }

    }

    // Return the integral.
    return integral / 64.0;
}

// The correct way to integrate a field over the mesh. The code is a modified
// version taken from the L2 projection mass matrix assembly, and uses Jacobians.
double integrateFieldDTKNew( DataTransferKit::UnitTest::ReferenceHexMesh &mesh,
                             DataTransferKit::Field &field, int integration_order )
{
    // Get the cells.
    int space_dim = 3;

    // Get the data.
    auto space            = mesh.functionSpace();
    auto set              = space->entitySet();
    auto shape_function   = space->shapeFunction();
    auto local_map        = space->localMap();
    auto integration_rule = space->integrationRule();
    auto cells            = space->entitySet()->entityIterator( space_dim );

    DataTransferKit::IntegrationPoint ip;
    ip.d_physical_coordinates.resize( space_dim );

    // Calculate the local integral.
    double integral  = 0.0;

    auto cells_begin = cells.begin();
    auto cells_end   = cells.end();
    for ( auto it = cells_begin; it != cells_end; ++it )
    {
        // Get the support ids of the entity.
        Teuchos::Array<DataTransferKit::SupportId> support_ids;
        shape_function->entitySupportIds( *it, support_ids );
        int cardinality = support_ids.size();

        // Get physical coordinates of support
        Teuchos::Array<Teuchos::Array<double>>  support_coordinates(cardinality);
        for (int ni = 0; ni < cardinality; ni++) {
            DataTransferKit::Entity support;
            set->getEntity(support_ids[ni], 0, support);

            support_coordinates[ni].resize(space_dim);
            local_map->centroid( support, support_coordinates[ni]() );
        }

        // Get the measure of the entity.
        double entity_measure = local_map->measure( *it );

        // Get the integration rule.
        Teuchos::Array<Teuchos::Array<double>> int_points;
        Teuchos::Array<double>                 int_weights;
        integration_rule->getIntegrationRule( *it, integration_order,
                                              int_points, int_weights );
        int num_ip = int_weights.size();


        // Create new integration points.
        for ( int p = 0; p < num_ip; ++p )
        {
            // Add owner data.
            ip.d_owner_measure      = entity_measure;
            ip.d_owner_support_ids  = support_ids;
            ip.d_integration_weight = int_weights[p];

            // Evaluate the shape function and gradient.
            Teuchos::Array<double>                  shape_evals;
            Teuchos::Array<Teuchos::Array<double>>  shape_grads;
            shape_function->evaluateValue( *it, int_points[p](),
                                           shape_evals );
            shape_function->evaluateGradient( *it, int_points[p](),
                                           shape_grads );
            ip.d_owner_shape_evals = shape_evals;

            // Compute Jacobian
            double J[space_dim][space_dim];
            for (int i = 0; i < space_dim; i++)
                for (int j = 0; j < space_dim; j++) {
                    J[i][j] = 0.0;
                    for (int ni = 0; ni < cardinality; ++ni)
                        J[i][j] += support_coordinates[ni][i] * shape_grads[ni][j];
                }

            // Compute derminant.
            double det = 0.0;
            if (space_dim == 2) {
                det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
            } else if (space_dim == 3) {
                det = J[0][0]*J[1][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] + J[0][2]*J[1][0]*J[2][1]
                        - J[0][2]*J[1][1]*J[2][0] - J[0][0]*J[1][2]*J[2][1] - J[2][2]*J[1][0]*J[0][1];
            }
            det *= 8; // for some reason

            // Map the integration point to the physical frame of the range
            // entity.
            local_map->mapToPhysicalFrame(
                *it, int_points[p](), ip.d_physical_coordinates() );

            // Update the integral.
            for ( int ni = 0; ni < cardinality; ++ni )
                integral += int_weights[p] * det * testFunction(ip.d_physical_coordinates());
        }

    }

    // Return the integral.
    return integral / 64.0;
}

// The correct way to integrate a field over the mesh. It uses Jacobians!
double integrateFieldIntrepid( DataTransferKit::UnitTest::ReferenceHexMesh &mesh,
                               DataTransferKit::Field &field, int integration_order )
{
    // Get the cells.
    int space_dim = 3;

    // Get the data.
    auto space            = mesh.functionSpace();
    auto set              = space->entitySet();
    auto shape_function   = space->shapeFunction();
    auto local_map        = space->localMap();
    auto integration_rule = space->integrationRule();
    auto cells            = space->entitySet()->entityIterator( space_dim );

    Teuchos::Array<DataTransferKit::SupportId>  support_ids;
    Teuchos::Array<Teuchos::Array<double>>      int_points;
    Teuchos::Array<double>                      int_weights;
    Teuchos::Array<double>                      jac_dets;;
    Teuchos::Array<Teuchos::Array<double>>      shape_evals;

    DataTransferKit::IntegrationPoint ip;
    ip.d_physical_coordinates.resize( space_dim );

    // Compute number of cells.
    int num_cells = 0;
    for ( auto cell = cells.begin(); cell != cells.end(); ++cell )
        num_cells++;

    // Populate Intrepid data.
    shards::CellTopology elem_topo = shards::getCellTopologyData<shards::Hexahedron<8>>();
    int num_nodes_per_elem = elem_topo.getNodeCount();

    Teuchos::Array<int> coords_dim(3);
    coords_dim[0] = num_cells;
    coords_dim[1] = num_nodes_per_elem;
    coords_dim[2] = space_dim;

    Intrepid::FieldContainer<double> elem_coords(num_cells, num_nodes_per_elem, space_dim);
    Intrepid::FieldContainer<double> dofs       (num_cells, num_nodes_per_elem);

    int node_ids[num_nodes_per_elem];
    Teuchos::Array<double> coords( space_dim );
    for (int i = 0; i < num_cells; i++)
    {
        compute_elem_node_ids( i, node_ids );

        for (int j = 0; j < num_nodes_per_elem; j++) {
            DataTransferKit::Entity node;
            set->getEntity(node_ids[j], 0, node);

            local_map->centroid( node, coords() );

            for (int k = 0; k < space_dim; k++) {
                elem_coords(i, j, k) = coords[k];
            }
            dofs(i,j) = testFunction(coords);

        }
    }

    DataTransferKit::IntrepidCell intrepid_cell(elem_topo, integration_order);
    intrepid_cell.allocateCellState( elem_coords );
    intrepid_cell.updateCellState();

    // Integrate.
    Intrepid::FieldContainer<double> integrals(num_cells);
    intrepid_cell.integrate(dofs, integrals);
    double integral = 0;
    for (int i = 0; i < num_cells; i++)
        integral += integrals(i);

    // Return the integral.
    return integral;
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( L2ProjectionOperator, integration )
{
    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

    DTK_REQUIRE(comm->getSize() == 1);

    // Set the global problem bounds.
    double x_min = 0.0;
    double y_min = 0.0;
    double z_min = 0.0;
    double x_max = 3.1;
    double y_max = 5.2;
    double z_max = 8.3;

    double perturb = 0.4;

    // Create a target mesh and field.
    int num_sx = NX;
    int num_sy = NY;
    int num_sz = NZ;

    DataTransferKit::UnitTest::ReferenceHexMesh mesh( comm,
        x_min, x_max, num_sx, y_min, y_max, num_sy, z_min, z_max, num_sz, perturb );

    auto field = mesh.nodalField( 1 );

    // Put some data on the nodal field.
    auto local_map = mesh.functionSpace()->localMap();
    auto nodes     = mesh.functionSpace()->entitySet()->entityIterator(0);

    auto nodes_begin = nodes.begin();
    auto nodes_end   = nodes.end();
    Teuchos::Array<double> coords( 3 );
    for ( nodes = nodes.begin(); nodes != nodes.end(); ++nodes ) {
        local_map->centroid( *nodes, coords() );

        field->writeFieldData( nodes->id(), 0, testFunction( coords() ) );
    }

    // Compare two integrals
    int integration_order = 3;
    double integralExact      = (x_max - x_min)*(y_max - y_min)*(z_max - z_min);
    double integralOld        = integrateFieldDTKOld  ( mesh, *field, integration_order );
    double integralNew        = integrateFieldDTKNew  ( mesh, *field, integration_order );
    double integralIntrepid   = integrateFieldIntrepid( mesh, *field, integration_order );

    std::cout << "perturbation = " << perturb << std::endl;
    std::cout << std::scientific << std::setprecision(3) << std::endl;
    std::cout << "integral (exact)    = " << integralExact << std::endl;
    std::cout << "integral (DTK old)  = " << integralOld << ", diff = "
            << std::abs(integralOld - integralExact) << std::endl;
    std::cout << "integral (DTK new)  = " << integralNew << ", diff = "
            << std::abs(integralNew - integralExact) << std::endl;
    std::cout << "integral (Intrepid) = " << integralIntrepid << ", diff = "
            << std::abs(integralIntrepid - integralExact)  << std::endl;
}

//---------------------------------------------------------------------------//
// end tstIntegral.cpp
//---------------------------------------------------------------------------//
