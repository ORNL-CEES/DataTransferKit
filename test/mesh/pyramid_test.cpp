#include <cstdlib>

#include <Teuchos_UnitTestHarness.hpp>

#include <Shards_BasicTopologies.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>

TEUCHOS_UNIT_TEST( Pyramid, pyramid_test )
{
    // Create a pyramid.
    shards::CellTopology topology(
	shards::getCellTopologyData<shards::Pyramid<5> >() );
    Intrepid::FieldContainer<double> pyramid_nodes( 1, 5, 3 );

    // Node 1.
    pyramid_nodes( 0, 0, 0 ) = 0.0;
    pyramid_nodes( 0, 0, 1 ) = 0.0;
    pyramid_nodes( 0, 0, 2 ) = 0.0;

    // Node 2.
    pyramid_nodes( 0, 1, 0 ) = 1.0;
    pyramid_nodes( 0, 1, 1 ) = 0.0;
    pyramid_nodes( 0, 1, 2 ) = 0.0;

    // Node 3.
    pyramid_nodes( 0, 2, 0 ) = 1.0;
    pyramid_nodes( 0, 2, 1 ) = 1.0;
    pyramid_nodes( 0, 2, 2 ) = 0.0;

    // Node 4.
    pyramid_nodes( 0, 3, 0 ) = 0.0;
    pyramid_nodes( 0, 3, 1 ) = 1.0;
    pyramid_nodes( 0, 3, 2 ) = 0.0;

    // Node 5.
    pyramid_nodes( 0, 4, 0 ) = 0.0;
    pyramid_nodes( 0, 4, 1 ) = 0.0;
    pyramid_nodes( 0, 4, 2 ) = 1.0;

    // Create random points to map to the reference frame. Some will be in the
    // reference frame and some outside.
    int num_points = 1000;
    Intrepid::FieldContainer<double> point( 1, 3 );
    Intrepid::FieldContainer<double> reference_point( 1, 3 );
    for ( int i = 0; i < num_points; ++i )
    {
	point( 0, 0 ) = 1.5 * (double) std::rand() / RAND_MAX - 0.25;
	point( 0, 1 ) = 1.5 * (double) std::rand() / RAND_MAX - 0.25;
	point( 0, 2 ) = 1.5 * (double) std::rand() / RAND_MAX - 0.25;

	Intrepid::CellTools<double>::mapToReferenceFrame(
	    reference_point, point, pyramid_nodes, topology, 0 );
    }
}
