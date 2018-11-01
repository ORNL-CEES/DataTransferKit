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

#ifndef DTK_LINEAR_BVH_FWD_HPP
#define DTK_LINEAR_BVH_FWD_HPP

#include <DTK_Box.hpp>
#include <DTK_DetailsBoundingVolumeHierarchyImpl.hpp>
#include <DTK_DetailsNode.hpp>

#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>

namespace DataTransferKit
{

template <typename DeviceType>
class BoundingVolumeHierarchy
{
  public:
    using device_type = DeviceType;
    using bounding_volume_type = Box;
    using size_type = typename Kokkos::View<Node *, DeviceType>::size_type;

    BoundingVolumeHierarchy() = default; // build an empty tree

    template <typename Primitives>
    BoundingVolumeHierarchy( Primitives const &primitives );

    KOKKOS_INLINE_FUNCTION
    size_type size() const { return _leaf_nodes.extent( 0 ); }

    KOKKOS_INLINE_FUNCTION
    bool empty() const { return size() == 0; }

    KOKKOS_INLINE_FUNCTION
    bounding_volume_type bounds() const
    {
        // NOTE should default constructor initialize to an invalid geometry?
        if ( empty() )
            return bounding_volume_type();
        return getBoundingVolume( getRoot() );
    }

    template <typename Predicates, typename... Args>
    inline void query( Predicates const &predicates, Args &&... args ) const
    {
        // FIXME lame placeholder for concept check
        static_assert( Kokkos::is_view<Predicates>::value, "must pass a view" );
        using Tag = typename Predicates::value_type::Tag;
        Details::BoundingVolumeHierarchyImpl<DeviceType>::queryDispatch(
            Tag{}, *this, predicates, std::forward<Args>( args )... );
    }

    // private:
  public:
    Kokkos::View<Node *, DeviceType, Kokkos::MemoryUnmanaged> getLeafNodes()
    {
        return _leaf_nodes;
    }

    /**
     * Return the root node of the BVH.
     */
    KOKKOS_INLINE_FUNCTION Node *getRoot()
    {
        if ( empty() )
            return nullptr;
        return ( size() > 1 ? _internal_nodes : _leaf_nodes ).data();
    }
    KOKKOS_INLINE_FUNCTION Node const *getRoot() const
    {
        if ( empty() )
            return nullptr;
        return ( size() > 1 ? _internal_nodes : _leaf_nodes ).data();
    }
    /**
     * Return the bounding volume of a node of the BVH.
     */
    KOKKOS_INLINE_FUNCTION Box &getBoundingVolume( Node *node )
    {
        return node->bounding_box;
    }
    KOKKOS_INLINE_FUNCTION Box const &
    getBoundingVolume( Node const *node ) const
    {
        return node->bounding_box;
    }
    /**
     * Return true if the node is a leaf.
     */
    KOKKOS_INLINE_FUNCTION bool isLeaf( Node const *node ) const
    {
        return ( node->children.first == nullptr );
    }
    /**
     * Return the index of a leaf node. If the passed node is not a leaf node,
     * this will return garbage.
     */
    // FIXME: debate on whether we only need getLeafPermutationIndex for leaf
    // nodes, and what to do (if we need to) for internal.
    // TODO: we can also think about adding a DTK_REQUIRE check here for leaf
    // nodes. It would like if (leaf->children.first == nullptr), as only one
    // child is nullptr here.
    KOKKOS_INLINE_FUNCTION int getLeafPermutationIndex( Node const *leaf ) const
    {
        static_assert( sizeof( size_t ) == sizeof( Node * ),
                       "Conversion is a bad idea if these sizes do not match" );
        return reinterpret_cast<size_t>( leaf->children.second );
    }

  private:
    friend struct Details::TreeTraversal<DeviceType>;

    Kokkos::View<Node *, DeviceType> _leaf_nodes;
    Kokkos::View<Node *, DeviceType> _internal_nodes;
};

} // namespace DataTransferKit

#endif
