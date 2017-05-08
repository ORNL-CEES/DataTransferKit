#ifndef DTK_DETAILS_TREE_TRAVERSAL_HPP
#define DTK_DETAILS_TREE_TRAVERSAL_HPP

#include <details/DTK_DetailsAlgorithms.hpp>
#include <details/DTK_Predicate.hpp>

#include <DTK_BVHQuery.hpp>
#include <DTK_LinearBVH.hpp> // BVH

#include <functional>
#include <list>
#include <queue>
#include <stack>
#include <utility> // std::pair, std::make_pair

namespace DataTransferKit
{

struct Node;

namespace Details
{

// priority queue helper for nearest neighbor search
using Value = std::pair<Node const *, double>;

struct CompareDistance
{
    bool operator()( Value const &lhs, Value const &rhs )
    {
        // larger distance means lower priority
        return lhs.second > rhs.second;
    }
};

using PriorityQueue =
    std::priority_queue<Value, std::vector<Value>, CompareDistance>;

// helper for k nearest neighbors search
// container to store candidates throughout the search
struct SortedList
{
    using Value = std::pair<int, double>;
    using Container = std::list<Value>;
    SortedList( int k )
        : _maxsize( k ){};
    template <typename... Args>
    void emplace( Args &&... args )
    {
        if ( _sorted_list.size() < _maxsize )
            _sorted_list.emplace_back( std::forward<Args>( args )... );
        else
        {
            Value val( std::forward<Args>( args )... );
            if ( !_compare( val, _sorted_list.back() ) )
                throw std::runtime_error(
                    "Error in DTK::Impl::SortedList::emplace()" );
            _sorted_list.back() = val;
        }
        _sorted_list.sort( _compare );
    }
    bool empty() const { return _sorted_list.empty(); }
    bool full() const { return _sorted_list.size() == _maxsize; }
    Container::size_type size() const { return _sorted_list.size(); }
    Value const &back() const { return _sorted_list.back(); }
    Container _sorted_list;
    // closer objects stored at the front of the list
    std::function<bool( Value const &a, Value const &b )> _compare =
        []( Value const &a, Value const &b ) { return a.second < b.second; };
    Container::size_type _maxsize;
};

// There are two (related) families of search: one using a spatial predicate and
// one using nearest neighbours query (see boost::geometry::queries
// documentation).
template <typename NO, typename Predicate>
KOKKOS_FUNCTION void
spatial_query( BVH<NO> const bvh, Predicate const &predicate, int *indices,
               unsigned int &n_indices, unsigned int n_max_indices )
{
    // Allocate traversal stack from thread-local memory, and push nullptr to
    // indicate that there are no postponed nodes.
    Node const *stack[64];
    Node const **stack_ptr = stack;
    *stack_ptr++ = nullptr;

    Node const *node = BVHQuery<NO>::getRoot( bvh );
    unsigned int pos = 0;
    do
    {
        if ( BVHQuery<NO>::isLeaf( node ) )
        {
            indices[pos] = BVHQuery<NO>::getIndex( bvh, node );
            ++pos;
#if HAVE_DTK_DBC
            if ( pos > n_max_indices )
                printf( "Increase the size of indices array\n" );
#else
            (void)n_max_indices;
#endif
        }
        else
        {
            for ( Node const *child :
                  {node->children.first, node->children.second} )
            {
                if ( predicate( child ) )
                {
                    *stack_ptr++ = child; // Push
#if HAVE_DTK_DBC
                    if ( stack_ptr - stack >= 64 )
                        printf( "Increase the size of the stack\n" );
#endif
                }
            }
        }

        node = *--stack_ptr; // Pop
    } while ( node != nullptr );
    n_indices = pos;
}

template <typename NO, typename Predicate>
unsigned int query_dispatch( BVH<NO> const bvh, Predicate const &pred,
                             Kokkos::View<int *, typename NO::device_type> out,
                             SpatialPredicateTag )
{
    unsigned int constexpr n_max_indices = 1000;
    int indices[n_max_indices];
    unsigned int n_indices = 0;
    spatial_query( bvh, pred, indices, n_indices, n_max_indices );
    out = Kokkos::View<int *, typename NO::device_type>( "dummy", n_indices );
    for ( unsigned int i = 0; i < n_indices; ++i )
    {
        out[i] = indices[i];
    }
    return n_indices;
}

// query k nearest neighbours
template <typename NO>
KOKKOS_FUNCTION void nearest_query( BVH<NO> const bvh, Point const &query_point,
                                    int k, int *indices,
                                    unsigned int &n_indices )
{
#ifndef KOKKOS_ENABLE_CUDA
    SortedList candidate_list( k );

    PriorityQueue queue;
    // priority does not matter for the root since the node will be
    // processed directly and removed from the priority queue we don't even
    // bother computing the distance to it
    Node const *node = BVHQuery<NO>::getRoot( bvh );
    double node_distance = 0.0;
    queue.emplace( node, node_distance );

    double cutoff = Kokkos::ArithTraits<double>::max();
    while ( !queue.empty() && node_distance < cutoff )
    {
        // get the node that is on top of the priority list (i.e. is the
        // closest to the query point)
        std::tie( node, node_distance ) = queue.top();
        queue.pop();
        if ( BVHQuery<NO>::isLeaf( node ) )
        {
            if ( node_distance < cutoff )
            {
                // add leaf node to the candidate list
                candidate_list.emplace( BVHQuery<NO>::getIndex( bvh, node ),
                                        node_distance );
                // update cutoff if k neighbors are in the list
                if ( candidate_list.full() )
                    std::tie( std::ignore, cutoff ) = candidate_list.back();
            }
        }
        else
        {
            // insert children of the node in the priority list
            for ( Node const *child :
                  {node->children.first, node->children.second} )
            {
                double child_distance =
                    distance( query_point, child->bounding_box );
                queue.emplace( child, child_distance );
            }
        }
    }

    n_indices = candidate_list.size();
    unsigned int pos = 0;
    for ( auto const &elem : candidate_list._sorted_list )
    {
        indices[pos] = elem.first;
        ++pos;
    }
#endif
}

template <typename NO, typename Predicate>
int query_dispatch( BVH<NO> const bvh, Predicate const &pred,
                    Kokkos::View<int *, typename NO::device_type> out,
                    NearestPredicateTag )
{
    unsigned int constexpr n_max_indices = 1000;
    int indices[n_max_indices];
    unsigned int n_indices = 0;
    nearest_query( bvh, pred._query_point, pred._k, indices, n_indices );
    int const n = n_indices;
    out = Kokkos::View<int *, typename NO::device_type>( "dummy", n );
    for ( unsigned int i = 0; i < n_indices; ++i )
        out[i] = indices[i];

    return n;
}

} // end namespace Details
} // end namespace DataTransferKit

#endif
