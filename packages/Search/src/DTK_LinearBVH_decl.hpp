/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_LINEAR_BVH_DECL_HPP
#define DTK_LINEAR_BVH_DECL_HPP

#include <DTK_Box.hpp>
#include <DTK_DetailsBoundingVolumeHierarchyImpl.hpp>
#include <DTK_DetailsNode.hpp>
#include <DTK_DetailsSortUtils.hpp>
#include <DTK_DetailsTreeConstruction_decl.hpp>

#include "DTK_LinearBVH_fwd.hpp"

#include <DTK_DetailsTreeConstruction_decl.hpp>

namespace DataTransferKit
{

template <typename DeviceType>
using BVH = BoundingVolumeHierarchy<DeviceType>;

template <typename DeviceType>
template <typename Primitives>
BoundingVolumeHierarchy<DeviceType>::BoundingVolumeHierarchy(
    Primitives const &primitives )
    : _leaf_nodes( Kokkos::ViewAllocateWithoutInitializing( "leaf_nodes" ),
                   primitives.extent( 0 ) )
    , _internal_nodes(
          Kokkos::ViewAllocateWithoutInitializing( "internal_nodes" ),
          primitives.extent( 0 ) > 0 ? primitives.extent( 0 ) - 1 : 0 )
{
    // FIXME lame placeholder for concept check
    static_assert( Kokkos::is_view<Primitives>::value, "must pass a view" );

    if ( empty() )
    {
        return;
    }

    if ( size() == 1 )
    {
        Kokkos::View<size_t *, DeviceType> permutation_indices( "permute", 1 );
        Details::TreeConstruction<DeviceType>::initializeLeafNodes(
            permutation_indices, primitives, _leaf_nodes );
        return;
    }

    // determine the bounding box of the scene
    Details::TreeConstruction<DeviceType>::calculateBoundingBoxOfTheScene(
        primitives, getBoundingVolume( getRoot() ) );

    // calculate morton code of all objects
    auto const n = primitives.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> morton_indices(
        Kokkos::ViewAllocateWithoutInitializing( "morton" ), n );
    Details::TreeConstruction<DeviceType>::assignMortonCodes(
        primitives, morton_indices, getBoundingVolume( getRoot() ) );

    // sort them along the Z-order space-filling curve
    auto permutation_indices = Details::sortObjects( morton_indices );

    Details::TreeConstruction<DeviceType>::initializeLeafNodes(
        permutation_indices, primitives, _leaf_nodes );

    // generate bounding volume hierarchy
    Details::TreeConstruction<DeviceType>::generateHierarchy(
        morton_indices, _leaf_nodes, _internal_nodes );

    // calculate bounding box for each internal node by walking the hierarchy
    // toward the root
    Details::TreeConstruction<DeviceType>::calculateBoundingBoxes( *this );
}

} // namespace DataTransferKit

#endif
