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
 * \file   tstInterpolationOperator.cpp
 * \brief  Interpolation Operotor unit tests.
 */
//---------------------------------------------------------------------------//

#include "DTK_InterpolationOperator.hpp"

#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_UnitTestHarness.hpp>

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( InterpolationOperator_BaseTopology, basic,
                                   SC, LO, GO, NO )
{
    // Test types.
    using InterpolationOperator =
        DataTransferKit::InterpolationOperator<SC, LO, GO, NO>;
    using scalar_type = SC;
    using local_ordinal_type = LO;
    using global_ordinal_type = GO;
    using node_type = NO;
    using device_type = typename NO::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    typedef Kokkos::Experimental::DynRankView<double, execution_space>
        DynRankView;

    shards::CellTopology cell_topology(
        shards::getCellTopologyData<shards::Hexahedron<8>>() );
    Teuchos::RCP<Intrepid2::Basis<execution_space>> basis = Teuchos::rcp(
        new Intrepid2::Basis_HGRAD_HEX_C1_FEM<execution_space>() );
    Intrepid2::EFunctionSpace function_space =
        Intrepid2::EFunctionSpace::FUNCTION_SPACE_HGRAD;

    // Hexahedron [-5,5]^3
    unsigned int const n_cells = 1;
    unsigned int n_nodes = 8;
    unsigned int space_dim = 3;
    DynRankView cell_nodes( "cell_nodes", n_cells, n_nodes, space_dim );
    cell_nodes( 0, 0, 0 ) = -5.;
    cell_nodes( 0, 0, 1 ) = 5.;
    cell_nodes( 0, 0, 2 ) = -5.;
    cell_nodes( 0, 1, 0 ) = 5.;
    cell_nodes( 0, 1, 1 ) = 5.;
    cell_nodes( 0, 1, 2 ) = -5.;
    cell_nodes( 0, 2, 0 ) = 5.;
    cell_nodes( 0, 2, 1 ) = -5.;
    cell_nodes( 0, 2, 2 ) = -5.;
    cell_nodes( 0, 3, 0 ) = -5.;
    cell_nodes( 0, 3, 1 ) = -5.;
    cell_nodes( 0, 3, 2 ) = -5.;
    cell_nodes( 0, 4, 0 ) = -5.;
    cell_nodes( 0, 4, 1 ) = 5.;
    cell_nodes( 0, 4, 2 ) = 5.;
    cell_nodes( 0, 5, 0 ) = 5.;
    cell_nodes( 0, 5, 1 ) = 5.;
    cell_nodes( 0, 5, 2 ) = 5.;
    cell_nodes( 0, 6, 0 ) = 5.;
    cell_nodes( 0, 6, 1 ) = -5.;
    cell_nodes( 0, 6, 2 ) = 5.;
    cell_nodes( 0, 7, 0 ) = -5.;
    cell_nodes( 0, 7, 1 ) = -5.;
    cell_nodes( 0, 7, 2 ) = 5.;

    // Create the interpolation operator
    InterpolationOperator interpolation_operator( basis, cell_topology,
                                                  cell_nodes, function_space );

    DynRankView coefficients( "coefficients", n_cells, n_nodes );
    coefficients( 0, 0 ) = 0.;
    coefficients( 0, 1 ) = 1.;
    coefficients( 0, 2 ) = 2.;
    coefficients( 0, 3 ) = 3.;
    coefficients( 0, 4 ) = 4.;
    coefficients( 0, 5 ) = 5.;
    coefficients( 0, 6 ) = 6.;
    coefficients( 0, 7 ) = 7.;

    DynRankView phys_points( "phys_points", 2, space_dim );
    phys_points( 0, 0 ) = 0.;
    phys_points( 0, 1 ) = 0.;
    phys_points( 0, 2 ) = 0.;
    phys_points( 1, 0 ) = 1.;
    phys_points( 1, 1 ) = 1.;
    phys_points( 1, 2 ) = 1.;

    unsigned int const n_points = 2;
    DynRankView values( "values", n_cells, n_points );

    // Compute the interpolation
    interpolation_operator.apply( values, coefficients, phys_points );

    for ( unsigned int i = 0; i < n_points; ++i )
        std::cout << values( 0, i ) << " ";
    std::cout << std::endl;
}

// Template instantiations

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( InterpolationOperator_BaseTopology,  \
                                          basic, SCALAR, LO, GO, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
