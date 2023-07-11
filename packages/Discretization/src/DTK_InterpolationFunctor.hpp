/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_INTERPOLATION_FUNCTOR_HPP
#define DTK_INTERPOLATION_FUNCTOR_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
namespace Functor
{
template <typename Scalar, typename BasisType, typename DeviceType>
class Interpolation
{
  public:
    Interpolation( unsigned int const dim,
                   Kokkos::View<Coordinate **, DeviceType> reference_points,
                   Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids,
                   Kokkos::View<Scalar **, DeviceType> dof_values,
                   Kokkos::View<Scalar **, DeviceType> output )
        : _dim( dim )
        , _n_basis( cell_dofs_ids.extent( 1 ) )
        , _n_fields( dof_values.extent( 1 ) )
        , _basis_values( "basis_values", output.extent( 0 ), _n_basis, dim )
        , _reference_points( reference_points )
        , _cell_dofs_ids( cell_dofs_ids )
        , _dof_values( dof_values )
        , _output( output )
    {
        DTK_REQUIRE( _output.extent( 1 ) == dof_values.extent( 1 ) );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( int const i ) const
    {
        auto ref_point = Kokkos::subview( _reference_points, i, Kokkos::ALL() );
        auto basis_values =
            Kokkos::subview( _basis_values, i, Kokkos::ALL(), Kokkos::ALL() );
        BasisType::getValues( basis_values, ref_point );

        for ( unsigned int j = 0; j < _n_basis; ++j )
            for ( unsigned int d = 0; d < _dim; ++d )
                for ( unsigned int k = 0; k < _n_fields; ++k )
                    _output( i, k ) += basis_values( j, d ) *
                                       _dof_values( _cell_dofs_ids( i, j ), k );
    }

  private:
    unsigned int const _dim;
    unsigned int const _n_basis;
    unsigned int const _n_fields;
    Kokkos::DynRankView<Coordinate, DeviceType> _basis_values;
    Kokkos::View<Coordinate **, DeviceType> _reference_points;
    Kokkos::View<LocalOrdinal **, DeviceType> _cell_dofs_ids;
    Kokkos::View<Scalar **, DeviceType> _dof_values;
    Kokkos::View<Scalar **, DeviceType> _output;
};

template <typename Scalar, typename BasisType, typename DeviceType>
class HgradInterpolation
{
  public:
    HgradInterpolation(
        Kokkos::View<Coordinate **, DeviceType> reference_points,
        Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids,
        Kokkos::View<Scalar **, DeviceType> dof_values,
        Kokkos::View<Scalar **, DeviceType> output )
        : _n_basis( cell_dofs_ids.extent( 1 ) )
        , _n_fields( dof_values.extent( 1 ) )
        , _basis_values( "basis_values", output.extent( 0 ), _n_basis )
        , _reference_points( reference_points )
        , _cell_dofs_ids( cell_dofs_ids )
        , _dof_values( dof_values )
        , _output( output )
    {
        DTK_REQUIRE( _output.extent( 1 ) == dof_values.extent( 1 ) );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( int const i ) const
    {
        auto ref_point = Kokkos::subview( _reference_points, i, Kokkos::ALL() );
        auto basis_values = Kokkos::subview( _basis_values, i, Kokkos::ALL() );
        BasisType::getValues( basis_values, ref_point );

        for ( unsigned int j = 0; j < _n_basis; ++j )
            for ( unsigned int k = 0; k < _n_fields; ++k )
                _output( i, k ) += basis_values( j ) *
                                   _dof_values( _cell_dofs_ids( i, j ), k );
    }

  private:
    unsigned int _n_basis;
    unsigned int _n_fields;
    // We cannot use Scalar because in Basis_HGRAD_PYR_C1_FEM there is a
    // check that basis_values and ref_point have the same type.
    Kokkos::View<Coordinate **, DeviceType> _basis_values;
    Kokkos::View<Coordinate **, DeviceType> _reference_points;
    Kokkos::View<LocalOrdinal **, DeviceType> _cell_dofs_ids;
    Kokkos::View<Scalar **, DeviceType> _dof_values;
    Kokkos::View<Scalar **, DeviceType> _output;
};
} // namespace Functor
} // namespace DataTransferKit

#endif
