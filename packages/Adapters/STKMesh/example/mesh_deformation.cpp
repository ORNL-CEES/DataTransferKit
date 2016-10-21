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
 * \file   mesh_deformation.cpp
 * \author Stuart Slattery
 * \brief  Mesh deformation example.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include "DTK_C_API.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_MultiVector.hpp>

#include <Intrepid_FieldContainer.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <Ioss_SubSystem.h>
#include <init/Ionit_Initializer.h>

//---------------------------------------------------------------------------//
// Set displacement field data on the boundary.
void setBoundaryDisplacement( stk::mesh::BulkData &bulk_data,
                              stk::mesh::Selector &fixed_selector,
                              stk::mesh::Selector &moving_selector,
                              stk::mesh::Selector &interior_selector,
                              Teuchos::Array<double> &coord_field,
                              Teuchos::Array<double> &displace_field, double ux,
                              double uy, int &num_boundary, int &num_interior )
{
    // Get the fixed boundary nodes.
    stk::mesh::BucketVector fixed_part_buckets =
        fixed_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> fixed_part_nodes;
    stk::mesh::get_selected_entities( fixed_selector, fixed_part_buckets,
                                      fixed_part_nodes );
    int fixed_num_part_nodes = fixed_part_nodes.size();

    // Get the moving boundary nodes.
    stk::mesh::BucketVector moving_part_buckets =
        moving_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> moving_part_nodes;
    stk::mesh::get_selected_entities( moving_selector, moving_part_buckets,
                                      moving_part_nodes );
    int moving_num_part_nodes = moving_part_nodes.size();

    // Get the number of boundary nodes.
    num_boundary = fixed_num_part_nodes + moving_num_part_nodes;

    // Get the interior nodes.
    stk::mesh::BucketVector interior_part_buckets =
        interior_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> interior_part_nodes;
    stk::mesh::get_selected_entities( interior_selector, interior_part_buckets,
                                      interior_part_nodes );
    num_interior = interior_part_nodes.size();

    // Allocate
    int total_num_nodes = num_boundary + num_interior;
    coord_field.resize( 3 * total_num_nodes );
    displace_field.resize( 3 * total_num_nodes );
    // NOTE: The mesh is actually 3D.

    // Get the coordinate field.
    const stk::mesh::FieldBase *coord_data_base =
        bulk_data.mesh_meta_data().coordinate_field();
    const stk::mesh::Field<double, stk::mesh::Cartesian3d> *coord_data =
        dynamic_cast<const stk::mesh::Field<double, stk::mesh::Cartesian3d> *>(
            coord_data_base );
    double *coords;

    // Create an offset.  f is an offset here. its used to index into a global
    // array. The coordinate array is broken up into chunks of the fixed
    // boundary, moving boundary, and interior node. The coordinates and
    // displacements for these exist in different STK mesh parts so we break
    // the coordinate and field assignment into loops over nodes in each part.
    int f = 0;

    // Get the fixed data.
    for ( int n = 0; n < fixed_num_part_nodes; ++n )
    {
        f = n;

        displace_field[f * 3] = 0.0;
        displace_field[f * 3 + 1] = 0.0;

        coords = stk::mesh::field_data( *coord_data, fixed_part_nodes[n] );
        coord_field[f * 3] = coords[0];
        coord_field[f * 3 + 1] = coords[1];
    }

    // Get the moving data.
    for ( int n = 0; n < moving_num_part_nodes; ++n )
    {
        f = n + fixed_num_part_nodes;

        displace_field[f * 3] = ux;
        displace_field[f * 3 + 1] = uy;

        coords = stk::mesh::field_data( *coord_data, moving_part_nodes[n] );
        coord_field[f * 3] = coords[0];
        coord_field[f * 3 + 1] = coords[1];
    }

    // Get the interior data.
    for ( int n = 0; n < num_interior; ++n )
    {
        f = n + num_boundary;

        displace_field[f * 3] = 0.0;
        displace_field[f * 3 + 1] = 0.0;

        coords = stk::mesh::field_data( *coord_data, interior_part_nodes[n] );
        coord_field[f * 3] = coords[0];
        coord_field[f * 3 + 1] = coords[1];
    }
}

//---------------------------------------------------------------------------//
// Deform the mesh.
void deformMesh( stk::mesh::BulkData &bulk_data,
                 stk::mesh::Selector &fixed_selector,
                 stk::mesh::Selector &moving_selector,
                 stk::mesh::Selector &interior_selector,
                 Teuchos::Array<double> &displace_field )
{
    // Get the fixed boundary nodes.
    stk::mesh::BucketVector fixed_part_buckets =
        fixed_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> fixed_part_nodes;
    stk::mesh::get_selected_entities( fixed_selector, fixed_part_buckets,
                                      fixed_part_nodes );
    int fixed_num_part_nodes = fixed_part_nodes.size();

    // Get the moving boundary nodes.
    stk::mesh::BucketVector moving_part_buckets =
        moving_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> moving_part_nodes;
    stk::mesh::get_selected_entities( moving_selector, moving_part_buckets,
                                      moving_part_nodes );
    int moving_num_part_nodes = moving_part_nodes.size();

    // Get the number of boundary nodes.
    int num_boundary_nodes = fixed_num_part_nodes + moving_num_part_nodes;

    // Get the interior nodes.
    stk::mesh::BucketVector interior_part_buckets =
        interior_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> interior_part_nodes;
    stk::mesh::get_selected_entities( interior_selector, interior_part_buckets,
                                      interior_part_nodes );
    int interior_num_part_nodes = interior_part_nodes.size();

    // Get the coordinate field.
    const stk::mesh::FieldBase *coord_field_base =
        bulk_data.mesh_meta_data().coordinate_field();
    const stk::mesh::Field<double, stk::mesh::Cartesian3d> *coord_field =
        dynamic_cast<const stk::mesh::Field<double, stk::mesh::Cartesian3d> *>(
            coord_field_base );
    double *coords;

    // Get the fixed boundary node coordinates.
    int f = 0;
    for ( int n = 0; n < fixed_num_part_nodes; ++n )
    {
        f = n;

        coords = stk::mesh::field_data( *coord_field, fixed_part_nodes[n] );
        coords[0] += displace_field[f * 3];
        coords[1] += displace_field[f * 3 + 1];
    }
    // NOTE: Recall mesh is 3D despite what the filename might say...

    // Get the moving boundary node coordinates.
    for ( int n = 0; n < moving_num_part_nodes; ++n )
    {
        f = n + fixed_num_part_nodes;

        coords = stk::mesh::field_data( *coord_field, moving_part_nodes[n] );
        coords[0] += displace_field[f * 3];
        coords[1] += displace_field[f * 3 + 1];
    }

    // Get the interior node coordinates.
    for ( int n = 0; n < interior_num_part_nodes; ++n )
    {
        f = n + num_boundary_nodes;

        coords = stk::mesh::field_data( *coord_field, interior_part_nodes[n] );
        coords[0] += displace_field[f * 3];
        coords[1] += displace_field[f * 3 + 1];
    }
}

//---------------------------------------------------------------------------//
// Example driver.
//
// This example takes a mesh (mesh_deform_2d.exo) and deforms it based on the
// movement of nodes on the interior boundary using spline interpolation. The
// exterior boundary nodes of the mesh are fixed. This technique is originally
// described in:
//
// A. Boer, M. vand er Schoot, and H. Bijl, "Mesh deformation based on radial
// basis function interpolation", Computers and Structures, vol 85,
// pp. 784-795, 2007.
//
// To execute the driver run:
//
// ./DataTransferKitSTKMeshAdapters_STKMeshDeformation.exe
//
//---------------------------------------------------------------------------//
int main( int argc, char *argv[] )
{
    // INITIALIZATION
    // --------------

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession( &argc, &argv );
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();

    // Set the mesh i/o data.
    std::string mesh_input_file = "mesh_deform_2d.exo";
    std::string mesh_output_file = "mesh_deform_2d_out.exo";
    std::string part_name = "domain";
    std::string fixed_bnd_name = "fixed_boundary";
    std::string moving_bnd_name = "moving_boundary";

    // Get the raw mpi communicator (basic typedef in STK).
    Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
        mpi_comm->getRawMpiComm();
    stk::ParallelMachine parallel_machine = ( *opaque_comm )();

    // MESH READ
    // ----------------

    // Load the mesh.
    stk::io::StkMeshIoBroker broker( parallel_machine );
    std::size_t input_index = broker.add_mesh_database(
        mesh_input_file, "exodus", stk::io::READ_MESH );
    broker.set_active_mesh( input_index );
    broker.create_input_mesh();

    // Get the domain part.
    stk::mesh::Part *part = broker.meta_data().get_part( part_name );
    stk::mesh::Selector global_selector( *part );

    // Create the source bulk data.
    broker.populate_bulk_data();
    Teuchos::RCP<stk::mesh::BulkData> bulk_data =
        Teuchos::rcpFromRef( broker.bulk_data() );

    // DISPLACEMENT SETUP
    // ------------------

    // Fixed boundary does not move.
    stk::mesh::Part *fixed_bnd_part =
        broker.meta_data().get_part( fixed_bnd_name );
    stk::mesh::Selector fixed_selector( *fixed_bnd_part );

    // Moving boundary moves.
    stk::mesh::Part *moving_bnd_part =
        broker.meta_data().get_part( moving_bnd_name );
    stk::mesh::Selector moving_selector( *moving_bnd_part );

    // Get the interior.
    stk::mesh::Selector interior_selector =
        global_selector - fixed_selector - moving_selector;

    // Define the displacement for the interior boundary. The exterior
    // boundary will not move.
    double ux = 1.0;
    double uy = 1.0;

    // Set the displacement.
    Teuchos::Array<double> coord_field;
    Teuchos::Array<double> displace_field;
    int num_boundary = 0;
    int num_interior = 0;
    setBoundaryDisplacement( *bulk_data, fixed_selector, moving_selector,
                             interior_selector, coord_field, displace_field, ux,
                             uy, num_boundary, num_interior );

    // DISPLACEMENT TRANSFER SETUP
    // ---------------------------

    // Use spline interpolation and a fairly large radius for a smooth
    // deformation of the mesh.
    std::string options =
        "{ \"Map Type\": \"Spline Interpolation\", \"RBF Radius\": \"4.0\" }";

    // Create the map between the boundary nodes and the interior nodes.
    DTK_Map *map = DTK_Map_create(
        parallel_machine, coord_field( 0, 3 * num_boundary ).getRawPtr(),
        num_boundary, DTK_INTERLEAVED,
        coord_field( 3 * num_boundary, 3 * num_interior ).getRawPtr(),
        num_interior, DTK_INTERLEAVED, 3, options.c_str() );
    // NOTE: The mesh is 3D!

    // DISPLACEMENT TRANSFER
    // ---------------------

    // Apply the map to transfer the displacement from the boundary nodes to
    // the interior nodes.
    DTK_Map_apply(
        map, displace_field( 0, 3 * num_boundary ).getRawPtr(), DTK_INTERLEAVED,
        displace_field( 3 * num_boundary, 3 * num_interior ).getRawPtr(),
        DTK_INTERLEAVED, 3, false );

    // cleanup
    DTK_Map_delete( map );

    // MESH DEFORMATION
    // ----------------

    // Apply the deformed coordinate field to the STK mesh
    deformMesh( *bulk_data, fixed_selector, moving_selector, interior_selector,
                displace_field );

    // DEFORMED MESH WRITE
    // -------------------

    // Write the STK mesh to file. Visualize to see the deformation.
    std::size_t output_index =
        broker.create_output_mesh( mesh_output_file, stk::io::WRITE_RESULTS );
    broker.begin_output_step( output_index, 0.0 );
    broker.end_output_step( output_index );
}

//---------------------------------------------------------------------------//
// end tstSTK_Mesh.cpp
//---------------------------------------------------------------------------//
