//---------------------------------------------------------------------------//
/*!
 * \file cxx_main.cpp
 * \author Stuart R. Slattery
 * \brief Geometry/geometry transfer example.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <algorithm> 
#include <cassert>

#include <DTK_VolumeSourceMap.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_FieldTools.hpp>
#include <DTK_CommTools.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Cylinder.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_Ptr.hpp>

//---------------------------------------------------------------------------//
// Source geometry creation.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Cylinder,int> >
createSourceGeometry( const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    int src_rank = comm->getRank();

    // Allocate space for 2 cylinders.
    Teuchos::ArrayRCP<DataTransferKit::Cylinder> cylinders( 2 );
    Teuchos::ArrayRCP<int> cylinder_ids( 2 );

    // Create 2 cylinders. These will be unique for each src proc.
    double length = 1.0;
    double radius = 0.25;

    cylinders[0] = DataTransferKit::Cylinder( length, radius, src_rank, 0, 0 );
    cylinder_ids[0] = src_rank;

    cylinders[1] = DataTransferKit::Cylinder( length, radius, src_rank+4, 0, 0 );
    cylinder_ids[1] = src_rank + 4;

    return Teuchos::rcp( 
	new DataTransferKit::GeometryManager<DataTransferKit::Cylinder,int>(
	    cylinders, cylinder_ids, comm, 3 ) );
}

//---------------------------------------------------------------------------//
// Target geometry creation.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Cylinder,int> >
createTargetGeometry( const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    int tgt_inv_rank = comm->getSize() - comm->getRank() - 1;

    // Allocate space for 2 cylinders.
    Teuchos::ArrayRCP<DataTransferKit::Cylinder> cylinders( 2 );
    Teuchos::ArrayRCP<int> cylinder_ids( 2 );

    // Create 2 cylinders. These will be unique for each tgt proc.
    double length = 1.0;
    double radius = 0.25;

    cylinders[0] = DataTransferKit::Cylinder( length, radius, 2*tgt_inv_rank, 0, 0 );
    cylinder_ids[0] = 2*tgt_inv_rank;

    cylinders[1] = DataTransferKit::Cylinder( length, radius, 2*tgt_inv_rank+1, 0, 0 );
    cylinder_ids[1] = 2*tgt_inv_rank + 1;

    return Teuchos::rcp( 
	new DataTransferKit::GeometryManager<DataTransferKit::Cylinder,int>(
	    cylinders, cylinder_ids, comm, 3 ) );
}

//---------------------------------------------------------------------------//
// Source function evaluator.
//---------------------------------------------------------------------------//
class SourceEvaluator: public DataTransferKit::FieldEvaluator<
    int, DataTransferKit::FieldContainer<double> >
{
  public:

    SourceEvaluator( const Teuchos::ArrayRCP<int>& geom_gids )
	: d_geom_gids( geom_gids )
    { 
	std::sort( geom_gids.begin(), geom_gids.end() );
    }

    ~SourceEvaluator()
    { /* ... */ }

    DataTransferKit::FieldContainer<double> evaluate( 
	const Teuchos::ArrayRCP<int>& gids,
	const Teuchos::ArrayRCP<double>& coords )
    {
	Teuchos::ArrayRCP<double> evaluated_data( gids.size() );
	for ( int n = 0; n < gids.size(); ++n )
	{
	    if ( std::binary_search( d_geom_gids.begin(), d_geom_gids.end(),
				     gids[n] ) )
	    {
		evaluated_data[n] = gids[n]+1;
	    }
	    else
	    {
		evaluated_data[n] = 0.0;
	    }
	}

	return DataTransferKit::FieldContainer<double>( evaluated_data, 1 );
    }

  private:

    Teuchos::ArrayRCP<int> d_geom_gids;
};

//---------------------------------------------------------------------------//
// Main function driver for the geometry/geometry transfer problem.
int main(int argc, char* argv[])
{
    // ---------------//
    // PARALLEL SETUP
    // ---------------//

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm_default = 
	Teuchos::DefaultComm<int>::getComm();
    int num_procs = comm_default->getSize();

    // Only 4 procs for this example. Full overlap between source and target.
    if( num_procs == 4 )
    {
	// ---------------//
	// SOURCE SETUP
	// ---------------//

	// Set required variables in the scope of the global communicator.
	Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Cylinder,int> >
	    src_geometry;
	Teuchos::RCP<
	    DataTransferKit::FieldEvaluator<int,DataTransferKit::FieldContainer<double> > >
	    src_evaluator;

	// Get source geometry.
	src_geometry = createSourceGeometry( comm_default );

	// Get the source field evaluator.
	src_evaluator = Teuchos::rcp( 
	    new SourceEvaluator( src_geometry->gids() ) );



	// ---------------//
	// TARGET SETUP
	// ---------------//

	// Set required variables in the scope of the global communicator.
	Teuchos::RCP<DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> > >
	    tgt_coords;
	Teuchos::ArrayRCP<double> tgt_array;
	Teuchos::RCP<
	    DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> > > 
	    tgt_space;

	Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Cylinder,int> >
	    tgt_geometry = createTargetGeometry( comm_default );

	// Extract the cylinder centroids as the blocked target coordinates.
	Teuchos::ArrayRCP<DataTransferKit::Cylinder> tgt_cylinders = 
	    tgt_geometry->geometry();
	int num_geom = tgt_cylinders.size();
	Teuchos::ArrayRCP<double> coords( num_geom*3 );
	for ( int i = 0; i < num_geom; ++i )
	{
	    Teuchos::Array<double> centroid = tgt_cylinders[i].centroid();
	    coords[i] = centroid[0];
	    coords[i+num_geom] = centroid[1];
	    coords[i+2*num_geom] = centroid[2];
	}
	Teuchos::RCP<DataTransferKit::FieldContainer<double> > coord_container = 
	    Teuchos::rcp( new DataTransferKit::FieldContainer<double>( coords, 3 ) );

	tgt_coords = Teuchos::rcp( 
	    new DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> >(
		coord_container, comm_default ) );

	// Get the target data space.
	tgt_array = Teuchos::ArrayRCP<double>( 2 );
	Teuchos::RCP<DataTransferKit::FieldContainer<double> > tgt_container = 
	    Teuchos::rcp( new DataTransferKit::FieldContainer<double>( tgt_array, 1 ) );
	tgt_space = Teuchos::rcp(
	    new DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> >(
		tgt_container, comm_default ) );
    


	// ---------------//
	// MAPPING SETUP
	// ---------------//

	// Create the mapping for the source-to-target transfer. The mapping will
	// occur over the union communicator in 3 dimensions. Keep track of missed
	// points as we expect to miss some.
	DataTransferKit::VolumeSourceMap<DataTransferKit::Cylinder,
					 int,
					 DataTransferKit::FieldContainer<double> >
	    src_to_tgt_map( comm_default, 3, true );

	// Setup the source-to-target map with the source geometry as the source
	// and the target coordinates as the target.
	src_to_tgt_map.setup( src_geometry, tgt_coords );



	// ---------------//
	// DATA TRANSFER
	// ---------------//
    
	// Apply the mapping to transfer the data.
	src_to_tgt_map.apply( src_evaluator, tgt_space );

	// Get the missed target centroids.
	Teuchos::ArrayView<int> missed = src_to_tgt_map.getMissedTargetPoints();

	// Check the resulting data transfer.
	int fail_count = 0;
	int tgt_inv_rank = comm_default->getSize() - comm_default->getRank() - 1;

	// First, check for missed points. There should be none.
	if ( missed.size() > 0 )
	{
	    ++fail_count;
	}

	// Next, check the data. 
	for ( int i = 0; i < 2; ++i )
	{
	    if ( tgt_array[i] != 1.0 + i + 2*tgt_inv_rank ) ++fail_count;
	}

	std::cout << std::endl;
	if ( fail_count == 0 )
	{
	    std::cout << "TEST PASSED: target proc " 
		      << comm_default->getRank() << std::endl;
	}
	else
	{
	    std::cout << "TEST FAILED " << fail_count << ": target proc " 
		      << comm_default->getRank() << std::endl;
	}
	std::cout << std::endl;

    }
    return 0;
}

//---------------------------------------------------------------------------//
// end cxx_main.hpp
//---------------------------------------------------------------------------//

