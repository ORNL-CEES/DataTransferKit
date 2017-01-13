#include <details/DTK_DetailsAlgorithms.hpp>
#include <details/DTK_DetailsTreeConstruction.hpp>

#include <algorithm>
#include <atomic>
#include <cassert> // TODO: probably want to use something else but fine for now

namespace DataTransferKit
{
namespace Details
{

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
unsigned int expandBits( unsigned int v )
{
    v = ( v * 0x00010001u ) & 0xFF0000FFu;
    v = ( v * 0x00000101u ) & 0x0F00F00Fu;
    v = ( v * 0x00000011u ) & 0xC30C30C3u;
    v = ( v * 0x00000005u ) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
unsigned int morton3D( double x, double y, double z )
{
    using std::min;
    using std::max;
    x = min( max( x * 1024.0, 0.0 ), 1023.0 );
    y = min( max( y * 1024.0, 0.0 ), 1023.0 );
    z = min( max( z * 1024.0, 0.0 ), 1023.0 );
    unsigned int xx = expandBits( (unsigned int)x );
    unsigned int yy = expandBits( (unsigned int)y );
    unsigned int zz = expandBits( (unsigned int)z );
    return xx * 4 + yy * 2 + zz;
}

// TODO: this is a mess
// we need a default impl
#define __clz( x ) __builtin_clz( x )
// default implementation if nothing else is available
// Taken from:
// http://stackoverflow.com/questions/23856596/counting-leading-zeros-in-a-32-bit-unsigned-integer-with-best-algorithm-in-c-pro
// WARNING: this implementation does **not** support __clz(0) (result should be
// 32 but this function returns 0)
int clz( uint32_t x )
{
    static const char debruijn32[32] = {
        0, 31, 9, 30, 3, 8,  13, 29, 2,  5,  7,  21, 12, 24, 28, 19,
        1, 10, 4, 14, 6, 22, 25, 20, 11, 15, 23, 26, 16, 27, 17, 18};
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x++;
    return debruijn32[x * 0x076be629 >> 27];
}

// TODO: use preprocessor directive to select an implementation
// it turns out NVDIA's implementation of int __clz(unsigned int x) is
// slightly different than GCC __builtin_clz
// this caused a bug in an early implementation of the function that compute
// the common prefixes between two keys (NB: when i == j)
int countLeadingZeros( unsigned int x )
{
#if defined __CUDACC__
    // intrinsic function that is only supported in device code
    // COMMENT: not sure how I am supposed to use it then...
    return __clz( x );

#elif defined __GNUC__
    // int __builtin_clz(unsigned int x) result is undefined if x is 0
    return x != 0 ? __builtin_clz( x ) : 32;
#else
    // similar problem with the default implementation
    return x != 0 ? clz( x ) : 32;
#endif
}

int commonPrefix( unsigned int const *k, int n, int i, int j )
{
    if ( j < 0 || j > n - 1 )
        return -1;
    // our construction algorithm relies on keys being unique so we handle
    // explicitly case of duplicate Morton codes by augmenting each key by a bit
    // representation of its index.
    if ( k[i] == k[j] )
    {
        // countLeadingZeros( k[i] ^ k[j] ) == 32
        return 32 + countLeadingZeros( i ^ j );
    }
    return countLeadingZeros( k[i] ^ k[j] );
}

// from "Thinking Parallel, Part III: Tree Construction on the GPU" by Karras
int findSplit( unsigned int *sorted_morton_codes, int first, int last )
{
    // Identical Morton codes => split the range in the middle.

    unsigned int first_code = sorted_morton_codes[first];
    unsigned int last_code = sorted_morton_codes[last];

    if ( first_code == last_code )
        return ( first + last ) >> 1;

    // Calculate the number of highest bits that are the same
    // for all objects, using the count-leading-zeros intrinsic.

    int common_prefix = __clz( first_code ^ last_code );

    // Use binary search to find where the next bit differs.
    // Specifically, we are looking for the highest object that
    // shares more than commonPrefix bits with the first one.

    int split = first; // initial guess
    int step = last - first;

    do
    {
        step = ( step + 1 ) >> 1;     // exponential decrease
        int new_split = split + step; // proposed new position

        if ( new_split < last )
        {
            unsigned int split_code = sorted_morton_codes[new_split];
            int split_prefix = __clz( first_code ^ split_code );
            if ( split_prefix > common_prefix )
                split = new_split; // accept proposal
        }
    } while ( step > 1 );

    return split;
}

// branchless sign function
// QUESTION: do we want to add it to DTK helpers?
inline int sgn( int x ) { return ( x > 0 ) - ( x < 0 ); }

Kokkos::pair<int, int> determineRange( unsigned int *sorted_morton_codes, int n,
                                       int i )
{
    using std::min;
    using std::max;
    // determine direction of the range (+1 or -1)
    int direction = sgn( commonPrefix( sorted_morton_codes, n, i, i + 1 ) -
                         commonPrefix( sorted_morton_codes, n, i, i - 1 ) );
    assert( direction == +1 || direction == -1 );
    // compute upper bound for the length of the range
    int max_step = 2;
    int common_prefix =
        commonPrefix( sorted_morton_codes, n, i, i - direction );
    // compute upper bound for the length of the range
    while ( commonPrefix( sorted_morton_codes, n, i,
                          i + direction * max_step ) > common_prefix )
    {
        max_step = max_step << 1;
    }
    // find the other end using binary search
    int split = 0;
    int step = max_step;
    do
    {
        step = step >> 1;
        if ( commonPrefix( sorted_morton_codes, n, i,
                           i + ( split + step ) * direction ) > common_prefix )
            split += step;
    } while ( step > 1 );
    int j = i + split * direction;
    return {min( i, j ), max( i, j )};
}

void calculateBoundingBoxOfTheScene( AABB const *bounding_boxes, int n,
                                     AABB &scene_bounding_box )
{
    // QUESTION: precondition on sceneBoundingBox?
    for ( int i = 0; i < n; ++i ) // parallel reduce
        expand( scene_bounding_box, bounding_boxes[i] );
}

// to assign the Morton code for a given object, we use the centroid point of
// its bounding box, and express it relative to the bounding box of the scene.
void assignMortonCodes( AABB const *bounding_boxes, unsigned int *morton_codes,
                        int n, AABB const &scene_bounding_box )
{
    Point xyz;
    double a, b;
    for ( int i = 0; i < n; ++i ) // parallel for
    {
        centroid( bounding_boxes[i], xyz );
        // scale coordinates with respect to bounding box of the scene
        for ( int d = 0; d < 3; ++d )
        {
            a = scene_bounding_box[2 * d + 0];
            b = scene_bounding_box[2 * d + 1];
            xyz[d] = ( xyz[d] - a ) / ( b - a );
        }
        morton_codes[i] = morton3D( xyz[0], xyz[1], xyz[2] );
    }
}

void sortObjects( unsigned int *morton_codes, int *object_ids, int n )
{
    using std::sort;
    // possibly use thrust::sort()
    // see https://thrust.github.io
    sort( object_ids, object_ids + n,
          [morton_codes]( int const &i, int const &j ) {
              return morton_codes[i] < morton_codes[j];
          } );
    // TODO: in-place permutation of mortonCodes rather than 2nd sort
    sort( morton_codes, morton_codes + n );
}

// from "Thinking Parallel, Part III: Tree Construction on the GPU" by Karras
Node *generateHierarchy( unsigned int *sorted_morton_codes, int n,
                         Node *leaf_nodes, Node *internal_nodes )
{
    // Construct internal nodes.

    for ( int i = 0; i < n - 1; ++i ) // in parallel
    {
        // Find out which range of objects the node corresponds to.
        // (This is where the magic happens!)

        auto range = determineRange( sorted_morton_codes, n, i );
        int first = range.first;
        int last = range.second;

        // Determine where to split the range.

        int split = findSplit( sorted_morton_codes, first, last );

        // Select childA.

        Node *childA;
        if ( split == first )
            childA = &leaf_nodes[split];
        else
            childA = &internal_nodes[split];

        // Select childB.

        Node *childB;
        if ( split + 1 == last )
            childB = &leaf_nodes[split + 1];
        else
            childB = &internal_nodes[split + 1];

        // Record parent-child relationships.

        internal_nodes[i].children.first = childA;
        internal_nodes[i].children.second = childB;
        childA->parent = &internal_nodes[i];
        childB->parent = &internal_nodes[i];
    }

    // Node 0 is the root.

    return &internal_nodes[0];
}

void calculateBoundingBoxes( Node const *leaf_nodes, Node *internal_nodes,
                             int n )
{
    // possibly use Kokkos::atomic_fetch_add() here
    std::vector<std::atomic_flag> atomic_flags( n - 1 );
    // flags are in an unspecified state on construction
    // their value cannot be copied/moved (constructor and assigment deleted)
    // so we have to loop over them and initialize them to the clear state
    for ( auto &flag : atomic_flags )
        flag.clear();

    Node *root = internal_nodes;
    for ( int i = 0; i < n; ++i ) // parallel for
    {
        Node *node = leaf_nodes[i].parent;
        while ( node != root )
        {
            if ( !atomic_flags[node - root].test_and_set() )
                break;
            for ( Node *child : {node->children.first, node->children.second} )
                expand( node->bounding_box, child->bounding_box );
            node = node->parent;
        }
        // NOTE: could stop at node != root and then just check that what we
        // computed earlier (bounding box of the scene) is indeed the union of
        // the two children.
    }
}

} // end namespace Details
} // end namespace DataTransferKit
