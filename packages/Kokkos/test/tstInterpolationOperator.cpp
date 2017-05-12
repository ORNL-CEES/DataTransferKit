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

#include "DTK_Basis.hpp"
#include "DTK_InterpolationOperator.hpp"
#include "DTK_Intrepid2Basis.hpp"

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Shards_CellTopologyManagedData.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_UnitTestHarness.hpp>

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( InterpolationOperator_BaseTopology_OneCell,
                                   basic, SC, LO, GO, NO )
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
    Teuchos::RCP<Intrepid2::Basis<execution_space>> basis_function =
        Teuchos::rcp(
            new Intrepid2::Basis_HGRAD_HEX_C1_FEM<execution_space>() );
    Intrepid2::EFunctionSpace function_space =
        Intrepid2::EFunctionSpace::FUNCTION_SPACE_HGRAD;
    Teuchos::RCP<DataTransferKit::Basis<SC, LO, GO, NO>> basis =
        Teuchos::rcp( new DataTransferKit::Intrepid2Basis<SC, LO, GO, NO>(
            basis_function, function_space, cell_topology ) );

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
    InterpolationOperator interpolation_operator( basis, cell_nodes );

    DynRankView coefficients( "coefficients", n_cells, n_nodes );
    coefficients( 0, 0 ) = 0.;
    coefficients( 0, 1 ) = 1.;
    coefficients( 0, 2 ) = 2.;
    coefficients( 0, 3 ) = 3.;
    coefficients( 0, 4 ) = 4.;
    coefficients( 0, 5 ) = 5.;
    coefficients( 0, 6 ) = 6.;
    coefficients( 0, 7 ) = 7.;

    unsigned int const n_points = 9;
    DynRankView phys_points( "phys_points", n_cells, n_points, space_dim );
    phys_points( 0, 0, 0 ) = 0.;
    phys_points( 0, 0, 1 ) = 0.;
    phys_points( 0, 0, 2 ) = 0.;
    phys_points( 0, 1, 0 ) = -5.;
    phys_points( 0, 1, 1 ) = 5.;
    phys_points( 0, 1, 2 ) = -5.;
    phys_points( 0, 2, 0 ) = 5.;
    phys_points( 0, 2, 1 ) = 5.;
    phys_points( 0, 2, 2 ) = -5.;
    phys_points( 0, 3, 0 ) = 5.;
    phys_points( 0, 3, 1 ) = -5.;
    phys_points( 0, 3, 2 ) = -5.;
    phys_points( 0, 4, 0 ) = -5.;
    phys_points( 0, 4, 1 ) = -5.;
    phys_points( 0, 4, 2 ) = -5.;
    phys_points( 0, 5, 0 ) = -5.;
    phys_points( 0, 5, 1 ) = 5.;
    phys_points( 0, 5, 2 ) = 5.;
    phys_points( 0, 6, 0 ) = 5.;
    phys_points( 0, 6, 1 ) = 5.;
    phys_points( 0, 6, 2 ) = 5.;
    phys_points( 0, 7, 0 ) = 5.;
    phys_points( 0, 7, 1 ) = -5.;
    phys_points( 0, 7, 2 ) = 5.;
    phys_points( 0, 8, 0 ) = -5.;
    phys_points( 0, 8, 1 ) = -5.;
    phys_points( 0, 8, 2 ) = 5.;

    // Compute the interpolation
    DynRankView values( "values", n_cells, n_points );
    interpolation_operator.apply( values, coefficients, phys_points );

    // Check the result
    std::vector<double> ref_values = {3.5, 0., 1., 2., 3., 4., 5., 6., 7.};
    TEST_COMPARE_ARRAYS( values, ref_values );
}

// Define our own topology and finite element
template <typename SC, typename LO, typename GO, typename NO>
class CustomBasis : public DataTransferKit::Basis<SC, LO, GO, NO>
{
  public:
    using scalar_type = SC;
    using local_ordinal_type = LO;
    using global_ordinal_type = GO;
    using node_type = NO;
    using device_type = typename NO::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    typedef Kokkos::Experimental::DynRankView<double, Kokkos::LayoutStride,
                                              execution_space>
        DynRankView;

    CustomBasis();

    void mapToReferenceFrame( DynRankView ref_points, DynRankView phys_points,
                              DynRankView cell_nodes ) override;

    // Basis functions on the reference cell [0,1] x [0,1] are:
    // xy, (1-x)*y, (1-x)*(1-y), and x*(1-y)
    void getValues( DynRankView ref_basis_values,
                    DynRankView const cell_ref_points ) override;

    unsigned int getCardinality() override;

    bool checkPointInclusion( DynRankView point ) override;

    Intrepid2::EFunctionSpace getEFunctionSpace() override;

  private:
    unsigned int const _n_dofs;
    Intrepid2::EFunctionSpace _function_space;
};

template <typename SC, typename LO, typename GO, typename NO>
CustomBasis<SC, LO, GO, NO>::CustomBasis()
    : _n_dofs( 4 )
    , _function_space( Intrepid2::EFunctionSpace::FUNCTION_SPACE_HGRAD )
{
}

// The geometry is a deformed unit square
template <typename SC, typename LO, typename GO, typename NO>
void CustomBasis<SC, LO, GO, NO>::mapToReferenceFrame(
    DynRankView ref_points, DynRankView phys_points,
    DynRankView /*cell_nodes*/ )
{
    unsigned int const n_points = phys_points.extent( 1 );
    for ( unsigned int i = 0; i < n_points; ++i )
    {
        double const x = phys_points( 0, i, 0 );
        ref_points( 0, i, 0 ) = x;
        ref_points( 0, i, 1 ) = phys_points( 0, i, 1 ) - x * x * x;
    }
}

template <typename SC, typename LO, typename GO, typename NO>
void CustomBasis<SC, LO, GO, NO>::getValues( DynRankView ref_basis_values,
                                             DynRankView const cell_ref_points )
{
    unsigned int const n_points = cell_ref_points.extent( 0 );
    for ( unsigned int i = 0; i < n_points; ++i )
    {
        double const x = cell_ref_points( i, 0 );
        double const y = cell_ref_points( i, 1 );
        ref_basis_values( 0, i ) = ( 1. - x ) * ( 1. - y );
        ref_basis_values( 1, i ) = x * ( 1. - y );
        ref_basis_values( 2, i ) = x * y;
        ref_basis_values( 3, i ) = ( 1. - x ) * y;
    }
}

template <typename SC, typename LO, typename GO, typename NO>
unsigned int CustomBasis<SC, LO, GO, NO>::getCardinality()
{
    return _n_dofs;
}

template <typename SC, typename LO, typename GO, typename NO>
bool CustomBasis<SC, LO, GO, NO>::checkPointInclusion( DynRankView point )
{
    return true;
}

template <typename SC, typename LO, typename GO, typename NO>
Intrepid2::EFunctionSpace CustomBasis<SC, LO, GO, NO>::getEFunctionSpace()
{
    return _function_space;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( InterpolationOperator_CustomTopology_OneCell,
                                   basic, SC, LO, GO, NO )
{
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

    Teuchos::RCP<DataTransferKit::Basis<SC, LO, GO, NO>> basis =
        Teuchos::rcp( new CustomBasis<SC, LO, GO, NO>() );

    unsigned int const n_cells = 1;
    unsigned int const n_nodes = 4;
    unsigned int const space_dim = 2;
    // Not used
    DynRankView cell_nodes( "cell_nodes", n_cells, n_nodes, space_dim );

    // Create the interpolation operator
    InterpolationOperator interpolation_operator( basis, cell_nodes );

    unsigned int const n_dofs = 4;
    DynRankView coefficients( "coefficients", n_cells, n_dofs );
    coefficients( 0, 0 ) = 0.;
    coefficients( 0, 1 ) = 1.;
    coefficients( 0, 2 ) = 2.;
    coefficients( 0, 3 ) = 3.;

    unsigned int const n_points = 6;
    DynRankView phys_points( "phys_points", n_cells, n_points, space_dim );
    phys_points( 0, 0, 0 ) = 0.;
    phys_points( 0, 0, 1 ) = 0.;
    phys_points( 0, 1, 0 ) = 1.;
    phys_points( 0, 1, 1 ) = 1.;
    phys_points( 0, 2, 0 ) = 1.;
    phys_points( 0, 2, 1 ) = 2.;
    phys_points( 0, 3, 0 ) = 0.;
    phys_points( 0, 3, 1 ) = 1.;
    phys_points( 0, 4, 0 ) = 0.5;
    phys_points( 0, 4, 1 ) = 0.125;
    phys_points( 0, 5, 0 ) = 0.5;
    phys_points( 0, 5, 1 ) = 1.125;

    // Compute the interpolation
    DynRankView values( "values", n_cells, n_points );
    interpolation_operator.apply( values, coefficients, phys_points );

    // Check the result
    std::vector<double> ref_values = {0., 1., 2., 3., 0.5, 2.5};
    TEST_COMPARE_ARRAYS( values, ref_values );
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
    InterpolationOperator_BaseTopology_MultipleCells, basic, SC, LO, GO, NO )
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
    Teuchos::RCP<Intrepid2::Basis<execution_space>> basis_function =
        Teuchos::rcp(
            new Intrepid2::Basis_HGRAD_HEX_C1_FEM<execution_space>() );
    Intrepid2::EFunctionSpace function_space =
        Intrepid2::EFunctionSpace::FUNCTION_SPACE_HGRAD;
    Teuchos::RCP<DataTransferKit::Basis<SC, LO, GO, NO>> basis =
        Teuchos::rcp( new DataTransferKit::Intrepid2Basis<SC, LO, GO, NO>(
            basis_function, function_space, cell_topology ) );

    unsigned int const n_cells = 2;
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
    cell_nodes( 1, 0, 0 ) = 5.;
    cell_nodes( 1, 0, 1 ) = 5.;
    cell_nodes( 1, 0, 2 ) = -5.;
    cell_nodes( 1, 1, 0 ) = 15.;
    cell_nodes( 1, 1, 1 ) = 5.;
    cell_nodes( 1, 1, 2 ) = -5.;
    cell_nodes( 1, 2, 0 ) = 15.;
    cell_nodes( 1, 2, 1 ) = -5.;
    cell_nodes( 1, 2, 2 ) = -5.;
    cell_nodes( 1, 3, 0 ) = 5.;
    cell_nodes( 1, 3, 1 ) = -5.;
    cell_nodes( 1, 3, 2 ) = -5.;
    cell_nodes( 1, 4, 0 ) = 5.;
    cell_nodes( 1, 4, 1 ) = 5.;
    cell_nodes( 1, 4, 2 ) = 5.;
    cell_nodes( 1, 5, 0 ) = 15.;
    cell_nodes( 1, 5, 1 ) = 5.;
    cell_nodes( 1, 5, 2 ) = 5.;
    cell_nodes( 1, 6, 0 ) = 15.;
    cell_nodes( 1, 6, 1 ) = -5.;
    cell_nodes( 1, 6, 2 ) = 5.;
    cell_nodes( 1, 7, 0 ) = 5.;
    cell_nodes( 1, 7, 1 ) = -5.;
    cell_nodes( 1, 7, 2 ) = 5.;

    // Create the interpolation operator
    InterpolationOperator interpolation_operator( basis, cell_nodes );

    DynRankView coefficients( "coefficients", n_cells, n_nodes );
    coefficients( 0, 0 ) = 0.;
    coefficients( 0, 1 ) = 1.;
    coefficients( 0, 2 ) = 2.;
    coefficients( 0, 3 ) = 3.;
    coefficients( 0, 4 ) = 4.;
    coefficients( 0, 5 ) = 5.;
    coefficients( 0, 6 ) = 6.;
    coefficients( 0, 7 ) = 7.;
    coefficients( 1, 0 ) = 0.;
    coefficients( 1, 1 ) = 1.;
    coefficients( 1, 2 ) = 2.;
    coefficients( 1, 3 ) = 3.;
    coefficients( 1, 4 ) = 4.;
    coefficients( 1, 5 ) = 5.;
    coefficients( 1, 6 ) = 6.;
    coefficients( 1, 7 ) = 7.;

    unsigned int const n_points = 1;
    DynRankView phys_points( "phys_points", n_cells, n_points, space_dim );
    phys_points( 0, 0, 0 ) = -5.;
    phys_points( 0, 0, 1 ) = 5.;
    phys_points( 0, 0, 2 ) = -5.;
    phys_points( 1, 0, 0 ) = 10.;
    phys_points( 1, 0, 1 ) = 0.;
    phys_points( 1, 0, 2 ) = 0.;

    // Compute the interpolation
    DynRankView values( "values", n_cells, n_points );
    interpolation_operator.apply( values, coefficients, phys_points );

    // Check the result
    std::vector<double> ref_values = {0., 3.5};
    TEST_COMPARE_ARRAYS( values, ref_values );
}

// Template instantiations

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                                      \
        InterpolationOperator_BaseTopology_OneCell, basic, SCALAR, LO, GO,     \
        NODE )                                                                 \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                                      \
        InterpolationOperator_CustomTopology_OneCell, basic, SCALAR, LO, GO,   \
        NODE )                                                                 \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(                                      \
        InterpolationOperator_BaseTopology_MultipleCells, basic, SCALAR, LO,   \
        GO, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
