#include <details/DTK_DetailsAlgorithms.hpp>
#include <details/DTK_DetailsTreeTraversal.hpp>

#include <Kokkos_ArithTraits.hpp>

namespace DataTransferKit
{
namespace Details
{

bool checkOverlap( AABB const &a, AABB const &b ) { return overlaps( a, b ); }

void traverseRecursive( CollisionList &list, BVH const &bvh,
                        const AABB &queryAABB, int queryObjectIdx,
                        Node const *node )
{
    // Bounding box overlaps the query => process node.
    if ( checkOverlap( node->bounding_box, queryAABB ) )
    {
        // Leaf node => report collision.
        if ( bvh.isLeaf( node ) )
            list.add( queryObjectIdx, bvh.getIndex( node ) );

        // Internal node => recurse to children.
        else
        {
            Node const *childL = node->children.first;
            Node const *childR = node->children.second;
            traverseRecursive( list, bvh, queryAABB, queryObjectIdx, childL );
            traverseRecursive( list, bvh, queryAABB, queryObjectIdx, childR );
        }
    }
}

void traverseIterative( CollisionList &list, BVH const &bvh,
                        AABB const &queryAABB, int queryObjectIdx )
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.
    Node const *stack[64];
    Node const **stackPtr = stack;
    *stackPtr++ = NULL; // push

    // Traverse nodes starting from the root.
    Node const *node = bvh.getRoot();
    do
    {
        // Check each child node for overlap.
        Node const *childL = node->children.first;
        Node const *childR = node->children.second;
        bool overlapL = ( checkOverlap( queryAABB, childL->bounding_box ) );
        bool overlapR = ( checkOverlap( queryAABB, childR->bounding_box ) );

        // Query overlaps a leaf node => report collision.
        if ( overlapL && bvh.isLeaf( childL ) )
            list.add( queryObjectIdx, bvh.getIndex( childL ) );

        if ( overlapR && bvh.isLeaf( childR ) )
            list.add( queryObjectIdx, bvh.getIndex( childR ) );

        // Query overlaps an internal node => traverse.
        bool traverseL = ( overlapL && !bvh.isLeaf( childL ) );
        bool traverseR = ( overlapR && !bvh.isLeaf( childR ) );

        if ( !traverseL && !traverseR )
            node = *--stackPtr; // pop
        else
        {
            node = ( traverseL ) ? childL : childR;
            if ( traverseL && traverseR )
                *stackPtr++ = childR; // push
        }
    } while ( node != NULL );
}

std::list<std::pair<int, double>>
within( BVH const &bvh, Point const &query_point, double radius )
{
    std::list<std::pair<int, double>> ret;

    Stack stack;

    Node const *node = bvh.getRoot();
    double node_distance = 0.0;
    stack.emplace( node, node_distance );
    while ( !stack.empty() )
    {
        std::tie( node, node_distance ) = stack.top();
        stack.pop();
        if ( bvh.isLeaf( node ) )
        {
            ret.emplace_back( bvh.getIndex( node ), node_distance );
        }
        else
        {
            for ( Node const *child :
                  {node->children.first, node->children.second} )
            {
                double child_distance =
                    distance( query_point, child->bounding_box );
                if ( child_distance <= radius )
                    stack.emplace( child, child_distance );
            }
        }
    }
    ret.sort(
        []( std::pair<int, double> const &a, std::pair<int, double> const &b ) {
            return a.second > b.second;
        } );
    return ret;
}

std::list<std::pair<int, double>> nearest( BVH const &bvh,
                                           Point const &query_point, int k )
{
    SortedList candidate_list( k );

    PriorityQueue queue;
    // priority does not matter for the root since the node will be processed
    // directly and removed from the priority queue
    // we don't even bother computing the distance to it
    Node const *node = bvh.getRoot();
    double node_distance = 0.0;
    queue.emplace( node, node_distance );

    double cutoff = Kokkos::ArithTraits<double>::max();
    while ( !queue.empty() && node_distance < cutoff )
    {
        // get the node that is on top of the priority list (i.e. is the closest
        // to the query point)
        std::tie( node, node_distance ) = queue.top();
        queue.pop();
        if ( bvh.isLeaf( node ) )
        {
            if ( node_distance < cutoff )
            {
                // add leaf node to the candidate list
                candidate_list.emplace( bvh.getIndex( node ), node_distance );
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
    return candidate_list._sorted_list;
}

} // end namespace Details
} // end namespace DataTransferKit
