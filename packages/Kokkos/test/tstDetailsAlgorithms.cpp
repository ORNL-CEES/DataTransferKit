#include <details/DTK_DetailsAlgorithms.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace dtk = DataTransferKit::Details;

TEUCHOS_UNIT_TEST( DetailsAlgorithms, distance )
{
    TEST_EQUALITY( dtk::distance( DataTransferKit::Point( {1.0, 2.0, 3.0} ),
                                  DataTransferKit::Point( {1.0, 1.0, 1.0} ) ),
                   std::sqrt( 5.0 ) );

    // box is unit cube
    DataTransferKit::BBox box( {0.0, 1.0, 0.0, 1.0, 0.0, 1.0} );
    // distance is zero if the point is inside the box
    TEST_EQUALITY(
        dtk::distance( DataTransferKit::Point( {0.5, 0.5, 0.5} ), box ), 0.0 );
    // or anywhere on the boundary
    TEST_EQUALITY(
        dtk::distance( DataTransferKit::Point( {0.0, 0.0, 0.5} ), box ), 0.0 );
    // normal projection onto center of one face
    TEST_EQUALITY(
        dtk::distance( DataTransferKit::Point( {2.0, 0.5, 0.5} ), box ), 1.0 );
    // projection onto edge
    TEST_EQUALITY(
        dtk::distance( DataTransferKit::Point( {2.0, 0.75, -1.0} ), box ),
        std::sqrt( 2.0 ) );
    // projection onto corner node
    TEST_EQUALITY(
        dtk::distance( DataTransferKit::Point( {-1.0, 2.0, 2.0} ), box ),
        std::sqrt( 3.0 ) );
}

TEUCHOS_UNIT_TEST( DetailsAlgorithms, overlaps )
{
    DataTransferKit::BBox box;
    // uninitialized box does not even overlap with itself
    TEST_ASSERT( !dtk::overlaps( box, box ) );
    // box with zero extent does
    box = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    TEST_ASSERT( dtk::overlaps( box, box ) );
    TEST_ASSERT( !dtk::overlaps( box, DataTransferKit::BBox() ) );
    // overlap with box that contains it
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0} ) ) );
    // does not overlap with some other box
    TEST_ASSERT( !dtk::overlaps(
        box, DataTransferKit::BBox( {1.0, 2.0, 1.0, 2.0, 1.0, 2.0} ) ) );
    // overlap when only touches another
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {0.0, 1.0, 0.0, 1.0, 0.0, 1.0} ) ) );
    // unit cube
    box = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    TEST_ASSERT( dtk::overlaps( box, box ) );
    TEST_ASSERT( !dtk::overlaps( box, DataTransferKit::BBox() ) );
    // smaller box inside
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {0.25, 0.75, 0.25, 0.75, 0.25, 0.75} ) ) );
    // bigger box that contains it
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {-1.0, 2.0, -1.0, 2.0, -1.0, 2.0} ) ) );
    // couple boxes that do overlap
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {0.5, 1.5, 0.5, 1.5, 0.5, 1.5} ) ) );
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {-0.5, 0.5, -0.5, 0.5, -0.5, 0.5} ) ) );
    // couple boxes that do not
    TEST_ASSERT( !dtk::overlaps(
        box, DataTransferKit::BBox( {-2.0, -1.0, -2.0, -1.0, -2.0, -1.0} ) ) );
    TEST_ASSERT( !dtk::overlaps(
        box, DataTransferKit::BBox( {0.0, 1.0, 0.0, 1.0, 2.0, 3.0} ) ) );
    // boxes overlap if faces touch
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {1.0, 2.0, 0.0, 1.0, 0.0, 1.0} ) ) );
    TEST_ASSERT( dtk::overlaps(
        box, DataTransferKit::BBox( {-0.5, 0.5, -0.5, 0.0, -0.5, 0.5} ) ) );
}

TEUCHOS_UNIT_TEST( DetailsAlgorithms, expand )
{
    // convenience utility to compare boxes
    auto checkBoxesEquality = []( DataTransferKit::BBox const &a,
                                  DataTransferKit::BBox const &b ) {
        for ( int i = 0; i < 6; ++i )
            if ( a[i] != b[i] )
                return false;
        return true;
    };
    // check that the utility does its job properly
    TEST_ASSERT( checkBoxesEquality(
        DataTransferKit::BBox( {0.0, 1.0, 0.0, 1.0, 0.0, 1.0} ),
        DataTransferKit::BBox( {0.0, 1.0, 0.0, 1.0, 0.0, 1.0} ) ) );
    TEST_ASSERT( !checkBoxesEquality(
        DataTransferKit::BBox( {0.0, 1.0, 0.0, 1.0, 0.0, 1.0} ),
        DataTransferKit::BBox( {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0} ) ) );

    DataTransferKit::BBox box;

    // expand box with points
    dtk::expand( box, DataTransferKit::Point( {0.0, 0.0, 0.0} ) );
    TEST_ASSERT( checkBoxesEquality(
        box, DataTransferKit::BBox( {0.0, 0.0, 0.0, 0.0, 0.0, 0.0} ) ) );
    dtk::expand( box, DataTransferKit::Point( {1.0, 1.0, 1.0} ) );
    TEST_ASSERT( checkBoxesEquality(
        box, DataTransferKit::BBox( {0.0, 1.0, 0.0, 1.0, 0.0, 1.0} ) ) );
    dtk::expand( box, DataTransferKit::Point( {0.25, 0.75, 0.25} ) );
    TEST_ASSERT( checkBoxesEquality(
        box, DataTransferKit::BBox( {0.0, 1.0, 0.0, 1.0, 0.0, 1.0} ) ) );
    dtk::expand( box, DataTransferKit::Point( {-1.0, -1.0, -1.0} ) );
    TEST_ASSERT( checkBoxesEquality(
        box, DataTransferKit::BBox( {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0} ) ) );

    // expand box with boxes
    dtk::expand(
        box, DataTransferKit::BBox( {0.25, 0.75, 0.25, 0.75, 0.25, 0.75} ) );
    TEST_ASSERT( checkBoxesEquality(
        box, DataTransferKit::BBox( {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0} ) ) );
    dtk::expand(
        box, DataTransferKit::BBox( {10.0, 11.0, 10.0, 11.0, 10.0, 11.0} ) );
    TEST_ASSERT( checkBoxesEquality(
        box, DataTransferKit::BBox( {-1.0, 11.0, -1.0, 11.0, -1.0, 11.0} ) ) );
}

TEUCHOS_UNIT_TEST( DetailsAlgorithms, centroid )
{
    DataTransferKit::BBox box( {-10.0, 0.0, 0.0, 10.0, 10.0, 20.0} );
    DataTransferKit::Point centroid;
    dtk::centroid( box, centroid );
    TEST_EQUALITY( centroid[0], -5.0 );
    TEST_EQUALITY( centroid[1], 5.0 );
    TEST_EQUALITY( centroid[2], 15.0 );
}
