#ifndef DTK_LINEAR_BVH_HPP
#define DTK_LINEAR_BVH_HPP

#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_Array.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_View.hpp>

#include <vector> // TODO: replace with Kokkos::View

namespace DataTransferKit
{

// Axis-Aligned Bounding Box
// This is just a thin wrapper around an array of size 2x spatial dimension with
// a default constructor to initialize properly an "empty" box.
struct AABB
{
    using ArrayType = Kokkos::Array<double, 6>;
    using SizeType = ArrayType::size_type;
    AABB() = default;
    AABB( ArrayType const &minmax )
        : _minmax( minmax )
    {
    }
    AABB &operator=( ArrayType const &minmax )
    {
        _minmax = minmax;
        return *this;
    }
    double &operator[]( SizeType i ) { return _minmax[i]; }
    double const &operator[]( SizeType i ) const { return _minmax[i]; }
    ArrayType _minmax = {{
        Kokkos::ArithTraits<double>::max(), -Kokkos::ArithTraits<double>::max(),
        Kokkos::ArithTraits<double>::max(), -Kokkos::ArithTraits<double>::max(),
        Kokkos::ArithTraits<double>::max(), -Kokkos::ArithTraits<double>::max(),
    }};
    friend std::ostream &operator<<( std::ostream &os, AABB const &aabb )
    {
        os << "{";
        for ( int d = 0; d < 3; ++d )
            os << " [" << aabb[2 * d + 0] << ", " << aabb[2 * d + 1] << "],";
        os << "}";
        return os;
    }
};
struct Node
{
    virtual ~Node() = default;
    Node *parent = nullptr;
    Kokkos::pair<Node *, Node *> children = {nullptr, nullptr};
    AABB bounding_box;
};
// Bounding Volume Hierarchy
struct BVH
{
    BVH( AABB const *bounding_boxes, int n );
    int size() const;
    AABB bounds() const;
    template <typename Predicate>
    int query( Predicate const &predicates, std::vector<int> &out ) const;
    bool isLeaf( Node const *node ) const;
    int getIndex( Node const *leaf ) const;
    Node const *getRoot() const;

  private:
    std::vector<Node> _leaf_nodes;
    std::vector<Node> _internal_nodes;
    // Array of indices that sort the boxes used to construct the hierarchy.
    // The leaf nodes are ordered so we need these to identify objects that meet
    // a predicate.
    std::vector<int> _indices;
};

} // end namespace DataTransferKit

#endif
