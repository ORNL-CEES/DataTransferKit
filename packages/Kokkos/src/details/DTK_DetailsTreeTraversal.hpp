#ifndef DTK_DETAILS_TREE_TRAVERSAL_HPP
#define DTK_DETAILS_TREE_TRAVERSAL_HPP

#include <details/DTK_DetailsAlgorithms.hpp> // overlap TODO:remove it

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

struct CollisionList
{
    void add( int i, int j ) { _ij.emplace_back( i, j ); }
    std::vector<std::pair<int, int>> _ij;
};
void traverseRecursive( CollisionList &list, BVH const &bvh,
                        AABB const &queryAABB, int queryObjectIdx,
                        Node const *node );
void traverseIterative( CollisionList &list, BVH const &bvh,
                        AABB const &queryAABB, int queryObjectIdx );
// TODO: get rid of this guy
bool checkOverlap( AABB const &a, AABB const &b );

// query k nearest neighbours
std::list<std::pair<int, double>>
nearest( BVH const &bvh, Point const &query_point, int k = 1 );

// radius search
std::list<std::pair<int, double>>
within( BVH const &bvh, Point const &query_point, double radius );

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
using Stack = std::stack<Value, std::vector<Value>>;

struct NearestPredicateTag
{
};
struct SpatialPredicateTag
{
};

template <typename Geometry>
struct Nearest
{
    using Tag = NearestPredicateTag;
    Nearest( Geometry const &geometry, int k )
        : _geometry( geometry )
        , _k( k )
    {
    }

    Geometry _geometry;
    int _k;
};

template <typename Geometry>
struct Within
{
    using Tag = SpatialPredicateTag;
    Within( Geometry const &geometry, double radius )
        : _geometry( geometry )
        , _radius( radius )
    {
    }

    Geometry _geometry;
    double _radius;
};

template <typename Geometry>
Nearest<Geometry> nearest( Geometry const &g, int k = 1 )
{
    return Nearest<Geometry>( g, k );
}

template <typename Geometry>
Within<Geometry> within( Geometry const &g, double r )
{
    return Within<Geometry>( g, r );
}

template <typename Predicate>
int query_dispatch( BVH const *bvh, Predicate const &pred,
                    std::vector<int> &out, NearestPredicateTag )
{
    auto nearest_neighbors = nearest( *bvh, pred._geometry, pred._k );
    int const n = nearest_neighbors.size();
    out.resize( n );
    int i = 0;
    for ( auto const &elem : nearest_neighbors )
    {
        out[i++] = elem.first;
    }
    if ( i != n )
        throw std::runtime_error( "ohhh no" );
    return n;
}

template <typename Predicate>
int query_dispatch( BVH const *bvh, Predicate const &pred,
                    std::vector<int> &out, SpatialPredicateTag )
{
    auto aaa = within( *bvh, pred._geometry, pred._radius );
    int const n = aaa.size();
    out.resize( n );
    int i = 0;
    for ( auto const &elem : aaa )
    {
        out[i++] = elem.first;
    }
    return n;
}

} // end namespace Details
} // end namespace DataTransferKit

#endif
