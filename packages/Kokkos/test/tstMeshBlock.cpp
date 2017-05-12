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
 * \file   tstMeshBlock.cpp
 * \brief  Mesh Block unit tests.
 */
//---------------------------------------------------------------------------//

#include <DTK_MeshBlock.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Kokkos_Core.hpp>

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MeshBlock, basic, SC, LO, GO, NO )
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
    using connectivity_view = Kokkos::View<local_ordinal_type **, device_type>;
    using const_global_id_view =
        Kokkos::View<const global_ordinal_type *, device_type>;
    using const_coordinate_view =
        Kokkos::View<const scalar_type **, device_type>;
    using const_connectivity_view =
        Kokkos::View<const local_ordinal_type **, device_type>;

    // Local numbering:
    // 6 -- 7 -- 8
    // |    |    |
    // 3 -- 4 -- 5
    // |    |    |
    // 0 -- 1 -- 2

    // Global numbering:
    // 2 -- 1 -- 0
    // |    |    |
    // 5 -- 4 -- 3
    // |    |    |
    // 8 -- 7 -- 6

    const unsigned int n_nodes = 9;
    global_id_view node_ids( "node_ids", n_nodes );
    for ( unsigned int i = 0; i < n_nodes; ++i )
        node_ids( i ) = n_nodes - ( i + 1 );

    const unsigned int n_cells = 4;
    connectivity_view connectivity( "connectiviy", n_cells, 4 );
    connectivity( 0, 0 ) = 0;
    connectivity( 0, 1 ) = 1;
    connectivity( 0, 2 ) = 4;
    connectivity( 0, 3 ) = 3;
    connectivity( 1, 0 ) = 1;
    connectivity( 1, 1 ) = 2;
    connectivity( 1, 2 ) = 5;
    connectivity( 1, 3 ) = 4;
    connectivity( 2, 0 ) = 3;
    connectivity( 2, 1 ) = 4;
    connectivity( 2, 2 ) = 7;
    connectivity( 2, 3 ) = 6;
    connectivity( 3, 0 ) = 4;
    connectivity( 3, 1 ) = 5;
    connectivity( 3, 2 ) = 8;
    connectivity( 3, 3 ) = 7;

    const unsigned int space_dim = 2;
    coordinate_view coordinates( "coordinates", n_nodes, space_dim );
    coordinates( 0, 0 ) = 0.;
    coordinates( 0, 1 ) = 0.;
    coordinates( 1, 0 ) = 1.;
    coordinates( 1, 1 ) = 0.;
    coordinates( 2, 0 ) = 2.;
    coordinates( 2, 1 ) = 0.;
    coordinates( 3, 0 ) = 0.;
    coordinates( 3, 1 ) = 1.;
    coordinates( 4, 0 ) = 1.;
    coordinates( 4, 1 ) = 1.;
    coordinates( 5, 0 ) = 2.;
    coordinates( 5, 1 ) = 1.;
    coordinates( 6, 0 ) = 0.;
    coordinates( 6, 1 ) = 2.;
    coordinates( 7, 0 ) = 1.;
    coordinates( 7, 1 ) = 2.;
    coordinates( 8, 0 ) = 2.;
    coordinates( 8, 1 ) = 2.;

    shards::CellTopology topology =
        shards::getCellTopologyData<shards::Quadrilateral<4>>();

    auto comm = Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();
    DataTransferKit::MeshBlock<SC, LO, GO, NO> mesh_block(
        comm, node_ids, connectivity, coordinates, topology );

    TEST_EQUALITY( mesh_block.spaceDim(), space_dim );
    TEST_EQUALITY( mesh_block.numLocalCells(), n_cells );
    TEST_EQUALITY( mesh_block.numGlobalCells(), comm_size * n_cells );
    TEST_EQUALITY( mesh_block.numLocalNodes(), n_nodes );
    TEST_EQUALITY( mesh_block.numGlobalNodes(), comm_size * n_nodes );

    const_global_id_view mesh_block_node_ids = mesh_block.nodeIds();
    for ( unsigned int i = 0; i < n_nodes; ++i )
        TEST_EQUALITY( mesh_block_node_ids( i ), node_ids( i ) );

    const_connectivity_view mesh_block_connectivity = mesh_block.connectivity();
    for ( unsigned int i = 0; i < n_cells; ++i )
        for ( unsigned int j = 0; j < 4; ++j )
            TEST_EQUALITY( mesh_block_connectivity( i, j ),
                           connectivity( i, j ) );

    const_coordinate_view mesh_block_coordinates = mesh_block.coordinates();
    for ( unsigned int i = 0; i < n_nodes; ++i )
        for ( unsigned int j = 0; j < space_dim; ++j )
            TEST_EQUALITY( mesh_block_coordinates( i, j ),
                           coordinates( i, j ) );

    TEST_EQUALITY( mesh_block.topology().getKey(), topology.getKey() );
}

//----------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MeshBlock, basic, SCALAR, LO, GO,    \
                                          NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
