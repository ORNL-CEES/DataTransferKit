#ifndef DTK_NODE_HPP
#define DTK_NODE_HPP

#include <DTK_Box.hpp>
#include <Kokkos_Pair.hpp>

namespace DataTransferKit
{
struct Node
{
    KOKKOS_INLINE_FUNCTION
    Node()
        : parent( nullptr )
        , children( {nullptr, nullptr} )
    {
    }

    Node *parent = nullptr;
    Kokkos::pair<Node *, Node *> children;
    BBox bounding_box;
};
}

#endif
