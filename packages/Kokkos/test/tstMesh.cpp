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
 * \file   tstMesh.cpp
 * \brief  Mesh unit tests.
 */
//---------------------------------------------------------------------------//

#include <DTK_Mesh.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Kokkos_Core.hpp>

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Mesh, basic, SC, LO, GO, NO )
{
    using scalar_type = SC;
    using local_ordinal_type = LO;
    using global_ordinal_type = GO;
    using node_type = NO;
    using device_type = typename NO::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    using global_id_view = Kokkos::View<global_ordinal_type *, device_type>;
    using coordinate_view = Kokkos::View<scalar_type **, device_type>;
    using connectivity_view = Kokkos::View<local_ordinal_type *, device_type>;
    using connectivity_view_2d =
        Kokkos::View<const local_ordinal_type **, device_type>;
    using const_global_id_view =
        Kokkos::View<const global_ordinal_type *, device_type>;
    using const_coordinate_view =
        Kokkos::View<const scalar_type **, device_type>;
    using const_connectivity_view =
        Kokkos::View<const local_ordinal_type *, device_type>;

    /* Local numbering:
         4 ---- 5
        /|      |\
       / |      | \
      /  |      |  \
     0 - 1 ---- 2 - 3 */

    /* Global numbering:
         1 ---- 0
        /|      |\
       / |      | \
      /  |      |  \
     5 - 4 ---- 3 - 2 */

    const unsigned int n_nodes = 6;
    global_id_view node_ids( "nodes_ids", n_nodes );
    for ( unsigned int i = 0; i < n_nodes; ++i )
        node_ids( i ) = n_nodes - ( i + 1 );

    connectivity_view connectivity( "connectivity", 10 );
    // First triangle
    connectivity( 0 ) = 0;
    connectivity( 1 ) = 1;
    connectivity( 2 ) = 4;
    // Quadrilateral
    connectivity( 3 ) = 1;
    connectivity( 4 ) = 2;
    connectivity( 5 ) = 5;
    connectivity( 6 ) = 4;
    // Second triangle
    connectivity( 7 ) = 2;
    connectivity( 8 ) = 2;
    connectivity( 9 ) = 5;

    const unsigned int space_dim = 2;
    coordinate_view coordinates( "coordinates", n_nodes, space_dim );
    coordinates( 0, 0 ) = 0.;
    coordinates( 0, 1 ) = 0.;
    coordinates( 1, 0 ) = 1.;
    coordinates( 1, 1 ) = 0.;
    coordinates( 2, 0 ) = 2.;
    coordinates( 2, 1 ) = 0.;
    coordinates( 3, 0 ) = 3.;
    coordinates( 3, 1 ) = 0.;
    coordinates( 4, 0 ) = 1.;
    coordinates( 4, 1 ) = 1.;
    coordinates( 5, 0 ) = 2.;
    coordinates( 5, 1 ) = 1.;

    std::vector<shards::CellTopology> topology;
    topology.push_back( shards::getCellTopologyData<shards::Triangle<3>>() );
    topology.push_back(
        shards::getCellTopologyData<shards::Quadrilateral<4>>() );
    topology.push_back( shards::getCellTopologyData<shards::Triangle<3>>() );

    auto comm = Teuchos::DefaultComm<int>::getComm();
    DataTransferKit::Mesh<SC, LO, GO, NO> mesh( comm, node_ids, connectivity,
                                                coordinates, topology );

    const Teuchos::RCP<std::vector<DataTransferKit::MeshBlock<SC, LO, GO, NO>>>
        mesh_blocks = mesh.meshBlocks();

    TEST_EQUALITY( ( *mesh_blocks )[0].topology().getKey(),
                   topology[0].getKey() );
    TEST_EQUALITY( ( *mesh_blocks )[1].topology().getKey(),
                   topology[1].getKey() );
    TEST_EQUALITY( ( *mesh_blocks )[0].numLocalCells(), 2 );
    TEST_EQUALITY( ( *mesh_blocks )[1].numLocalCells(), 1 );

    connectivity_view_2d connectivity_triangle =
        ( *mesh_blocks )[0].connectivity();
    TEST_EQUALITY( connectivity_triangle( 0, 0 ), connectivity( 0 ) );
    TEST_EQUALITY( connectivity_triangle( 0, 1 ), connectivity( 1 ) );
    TEST_EQUALITY( connectivity_triangle( 0, 2 ), connectivity( 2 ) );
    TEST_EQUALITY( connectivity_triangle( 1, 0 ), connectivity( 7 ) );
    TEST_EQUALITY( connectivity_triangle( 1, 1 ), connectivity( 8 ) );
    TEST_EQUALITY( connectivity_triangle( 1, 2 ), connectivity( 9 ) );

    connectivity_view_2d connectivity_quad = ( *mesh_blocks )[1].connectivity();
    TEST_EQUALITY( connectivity_quad( 0, 0 ), connectivity( 3 ) );
    TEST_EQUALITY( connectivity_quad( 0, 1 ), connectivity( 4 ) );
    TEST_EQUALITY( connectivity_quad( 0, 2 ), connectivity( 5 ) );
    TEST_EQUALITY( connectivity_quad( 0, 3 ), connectivity( 6 ) );
}

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Mesh, basic, SCALAR, LO, GO, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
