#include <details/DTK_DetailsAlgorithms.hpp>
#include <details/DTK_DetailsTreeConstruction.hpp>

#include <Kokkos_ArithTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <algorithm>
#include <bitset>
#include <sstream>
#include <vector>

namespace dtk = DataTransferKit::Details;

TEUCHOS_UNIT_TEST( DetailsBVH, morton_codes )
{
    std::vector<dtk::Point> points = {
        dtk::Point( {0.0, 0.0, 0.0} ),
        dtk::Point( {0.25, 0.75, 0.25} ),
        dtk::Point( {0.75, 0.25, 0.25} ),
        dtk::Point( {0.75, 0.75, 0.25} ),
        dtk::Point( {1.33, 2.33, 3.33} ),
        dtk::Point( {1.66, 2.66, 3.66} ),
        dtk::Point( {1024.0, 1024.0, 1024.0} ),
    };
    int const n = points.size();
    // lower left front corner corner of the octant the points fall in
    std::vector<std::array<unsigned int, 3>> anchors = {
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0},         {0, 0, 0},
        {1, 2, 3}, {1, 2, 3}, {1023, 1023, 1023}};
    auto fun = []( std::array<unsigned int, 3> const &anchor ) {
        unsigned int i = std::get<0>( anchor );
        unsigned int j = std::get<1>( anchor );
        unsigned int k = std::get<2>( anchor );
        return 4 * dtk::expandBits( i ) + 2 * dtk::expandBits( j ) +
               dtk::expandBits( k );
    };
    std::vector<unsigned int> ref( n,
                                   Kokkos::ArithTraits<unsigned int>::max() );
    for ( int i = 0; i < n; ++i )
        ref[i] = fun( anchors[i] );
    // using points rather than boxes for convenience here but still have to
    // build the axis-aligned bounding boxes around them
    std::vector<dtk::Box> boxes( n );
    for ( int i = 0; i < n; ++i )
        dtk::expand( boxes[i], points[i] );

    dtk::Box scene;
    dtk::calculateBoundingBoxOfTheScene( boxes.data(), n, scene );
    for ( int d = 0; d < 3; ++d )
    {
        TEST_EQUALITY( scene[2 * d + 0], 0.0 );
        TEST_EQUALITY( scene[2 * d + 1], 1024.0 );
    }

    std::vector<unsigned int> morton_codes(
        n, Kokkos::ArithTraits<unsigned int>::max() );
    dtk::assignMortonCodes( boxes.data(), morton_codes.data(), n, scene );
    for ( int i = 0; i < n; ++i )
        TEST_EQUALITY( morton_codes[i], ref[i] );
}

TEUCHOS_UNIT_TEST( DetailsBVH, indirect_sort )
{
    // need a functionality that sort objects based on their Morton code and
    // also returns the indices in the original configuration

    // dummy unsorted Morton codes and corresponding sorted indices as reference
    // solution
    std::vector<unsigned int> k = {2, 1, 4, 3};
    std::vector<int> ref = {1, 0, 3, 2};
    // distribute ids to unsorted objects
    int const n = k.size();
    std::vector<int> ids( n );
    std::iota( ids.begin(), ids.end(), 0 );
    // sort morton codes and object ids
    dtk::sortObjects( k.data(), ids.data(), n );
    // check that they are sorted
    TEST_ASSERT( std::is_sorted( k.begin(), k.end() ) );
    // check that ids are properly ordered
    for ( int i = 0; i < n; ++i )
        TEST_EQUALITY( ids[i], ref[i] );
}

TEUCHOS_UNIT_TEST( DetailsBVH, number_of_leading_zero_bits )
{
    TEST_EQUALITY( dtk::countLeadingZeros( 0 ), 32 );
    TEST_EQUALITY( dtk::countLeadingZeros( 1 ), 31 );
    TEST_EQUALITY( dtk::countLeadingZeros( 2 ), 30 );
    TEST_EQUALITY( dtk::countLeadingZeros( 3 ), 30 );
    TEST_EQUALITY( dtk::countLeadingZeros( 4 ), 29 );
    TEST_EQUALITY( dtk::countLeadingZeros( 5 ), 29 );
    TEST_EQUALITY( dtk::countLeadingZeros( 6 ), 29 );
    TEST_EQUALITY( dtk::countLeadingZeros( 7 ), 29 );
    TEST_EQUALITY( dtk::countLeadingZeros( 8 ), 28 );
    TEST_EQUALITY( dtk::countLeadingZeros( 9 ), 28 );
    // bitwise exclusive OR operator to compare bits
    TEST_EQUALITY( dtk::countLeadingZeros( 1 ^ 0 ), 31 );
    TEST_EQUALITY( dtk::countLeadingZeros( 2 ^ 0 ), 30 );
    TEST_EQUALITY( dtk::countLeadingZeros( 2 ^ 1 ), 30 );
    TEST_EQUALITY( dtk::countLeadingZeros( 3 ^ 0 ), 30 );
    TEST_EQUALITY( dtk::countLeadingZeros( 3 ^ 1 ), 30 );
    TEST_EQUALITY( dtk::countLeadingZeros( 3 ^ 2 ), 31 );
    TEST_EQUALITY( dtk::countLeadingZeros( 4 ^ 0 ), 29 );
    TEST_EQUALITY( dtk::countLeadingZeros( 4 ^ 1 ), 29 );
    TEST_EQUALITY( dtk::countLeadingZeros( 4 ^ 2 ), 29 );
    TEST_EQUALITY( dtk::countLeadingZeros( 4 ^ 3 ), 29 );
}

TEUCHOS_UNIT_TEST( DetailsBVH, common_prefix )
{
    // NOTE: Morton codes below are **not** unique
    std::vector<unsigned int> const fi = {
        0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144,
    };
    int const n = fi.size();
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 0, 0 ), 32 + 32 );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 0, 1 ), 31 );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 1, 0 ), 31 );
    // duplicate Morton codes
    TEST_EQUALITY( fi[1], 1 );
    TEST_EQUALITY( fi[1], fi[2] );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 1, 1 ), 64 );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 1, 2 ), 32 + 30 );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 2, 1 ), 62 );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 2, 2 ), 64 );
    // by definition \delta(i, j) = -1 when j \notin [0, n-1]
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 0, -1 ), -1 );
    TEST_EQUALITY( n, 13 );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 12, 12 ), 64 );
    TEST_EQUALITY( dtk::commonPrefix( fi.data(), n, 12, 13 ), -1 );
}

TEUCHOS_UNIT_TEST( DetailsBVH, example_tree_construction )
{
    // This is the example from the articles by Karras.
    // See
    // https://devblogs.nvidia.com/parallelforall/thinking-parallel-part-iii-tree-construction-gpu/
    std::vector<unsigned int> sorted_morton_codes;
    for ( std::string const &s : {
              "00001", "00010", "00100", "00101", "10011", "11000", "11001",
              "11110",
          } )
    {
        std::bitset<6> b( s );
        std::cout << b << "  " << b.to_ulong() << "\n";
        sorted_morton_codes.push_back( b.to_ulong() );
    }
    int const n = sorted_morton_codes.size();

    // reference solution for a recursive traversal from top to bottom
    // starting from root, visiting first the left child and then the right one
    std::ostringstream ref;
    ref << "I0"
        << "I3"
        << "I1"
        << "L0"
        << "L1"
        << "I2"
        << "L2"
        << "L3"
        << "I4"
        << "L4"
        << "I5"
        << "I6"
        << "L5"
        << "L6"
        << "L7";
    std::cout << "ref=" << ref.str() << "\n";

    // hierarchy generation
    std::vector<DataTransferKit::Node> leaf_nodes( n );
    std::vector<DataTransferKit::Node> internal_nodes( n - 1 );
    std::function<void( DataTransferKit::Node *, std::ostream & )>
        traverseRecursive;
    traverseRecursive = [&leaf_nodes, &internal_nodes, &traverseRecursive](
        DataTransferKit::Node *node, std::ostream &os ) {
        if ( std::any_of( leaf_nodes.begin(), leaf_nodes.end(),
                          [node]( DataTransferKit::Node const &leaf_node ) {
                              return std::addressof( leaf_node ) == node;
                          } ) )
        {
            os << "L" << node - leaf_nodes.data();
        }
        else
        {
            os << "I" << node - internal_nodes.data();
            for ( DataTransferKit::Node *child :
                  {node->children.first, node->children.second} )
                traverseRecursive( child, os );
        }
    };

    DataTransferKit::Details::generateHierarchy( sorted_morton_codes.data(), n,
                                                 leaf_nodes.data(),
                                                 internal_nodes.data() );

    DataTransferKit::Node *root = internal_nodes.data();
    TEST_ASSERT( root->parent == nullptr );

    std::ostringstream sol;
    traverseRecursive( root, sol );
    std::cout << "sol=" << sol.str() << "\n";

    TEST_EQUALITY( sol.str().compare( ref.str() ), 0 );
}
