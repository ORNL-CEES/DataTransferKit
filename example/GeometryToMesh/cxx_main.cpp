//---------------------------------------------------------------------------//
/*!
 * \file cxx_main.cpp
 * \author Stuart R. Slattery
 * \brief Geometry/mesh transfer example.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <algorithm> 
#include <cassert>

#include <DTK_VolumeSourceMap.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraitsFieldAdapter.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_FieldTools.hpp>
#include <DTK_CommTools.hpp>
#include <DTK_CommIndexer.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Box.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_Ptr.hpp>

//---------------------------------------------------------------------------//
// Source geometry creation.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Box,int> >
createSourceGeometry( const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    int src_rank = comm->getRank();

    // Allocate space for 2 boxes.
    Teuchos::ArrayRCP<DataTransferKit::Box> boxes( 2 );
    Teuchos::ArrayRCP<int> box_ids( 2 );

    // Create 2 boxes. 1 box will be shared with the other proc.
    boxes[0] = DataTransferKit::Box( src_rank, 0, 0,
				     src_rank+1, 1, 1 );
    box_ids[0] = src_rank;

    boxes[1] = DataTransferKit::Box( src_rank+1, 0, 0,
				     src_rank+2, 1, 1 );
    box_ids[1] = src_rank + 1;

    return Teuchos::rcp( 
	new DataTransferKit::GeometryManager<DataTransferKit::Box,int>(
	    boxes, box_ids, comm, 3 ) );
}

//---------------------------------------------------------------------------//
// Target mesh creation. Shifted in +X by 0.5 compared to the source so that
// we miss points on purpose.
//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::MeshContainer<int> >
createTargetMesh( const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    int tgt_rank = comm->getRank();
    
    // Create the vertex ids.
    int num_verts = 12;
    Teuchos::ArrayRCP<int> vertex_ids( num_verts );
    for ( int i = 0; i < num_verts; ++i ) 
	vertex_ids[i] = tgt_rank*num_verts + i;

    // Create the vertex coordinates. Blocked.
    int num_coords = num_verts*3;
    Teuchos::ArrayRCP<double> vertex_coords( num_coords );
    for ( int k = 0; k < 2; ++k ) {
	for ( int j = 0; j < 2; ++j ) {
	    for ( int i = 0; i < 3; ++i ) {
		int index = i + 3*j + 6*k;
		vertex_coords[index] = (i + tgt_rank)*1.0 + 0.5;
		vertex_coords[index+num_verts] = j*1.0;
		vertex_coords[index+2*num_verts] = k*1.0;
	    }
	}
    }

    // Create the element ids.
    int num_elements = 2;
    Teuchos::ArrayRCP<int> element_ids( num_elements );
    element_ids[0] = tgt_rank*2;
    element_ids[1] = tgt_rank*2 + 1;

    // Create the element connectivity from the vertex ids.
    int conn_size = 8;
    Teuchos::ArrayRCP<int> element_conn( num_elements*conn_size );
    element_conn[0]  = 0;   element_conn[1]  = 1;
    element_conn[2]  = 1;   element_conn[3]  = 2;
    element_conn[4]  = 4;   element_conn[5]  = 5;
    element_conn[6]  = 3;   element_conn[7]  = 4;
    element_conn[8]  = 6;   element_conn[9]  = 7;
    element_conn[10] = 7;   element_conn[11] = 8;
    element_conn[12] = 10;  element_conn[13] = 11;
    element_conn[14] = 9;   element_conn[15] = 10;

    // Create the permutation list.
    Teuchos::ArrayRCP<int> permutation( conn_size );
    for ( int i = 0; i < conn_size; ++i ) permutation[i] = i;

    // Create the mesh container.
    return Teuchos::rcp( new DataTransferKit::MeshContainer<int>( 
			     3,
			     vertex_ids,
			     vertex_coords,
			     DataTransferKit::DTK_HEXAHEDRON,
			     conn_size,
			     element_ids,
			     element_conn,
			     permutation ) );
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
// Main function driver for the geometry/mesh transfer problem.
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

    // Only 4 procs for this example.
    assert( num_procs == 4 );

    // Split the main communicator into 2 separate groups, one for the source
    // geometry and one for the target mesh.
    Teuchos::Array<int> sub_ranks_src(2), sub_ranks_tgt(2);
    for ( int i = 0; i < 2; ++i )
    {
	sub_ranks_src[i] = i;
	sub_ranks_tgt[i] = i+2;
    }

    // Generate the source and target communicators from the sub ranks.
    Teuchos::RCP<const Teuchos::Comm<int> > src_comm = 
	comm_default->createSubcommunicator( sub_ranks_src() );

    Teuchos::RCP<const Teuchos::Comm<int> > tgt_comm = 
	comm_default->createSubcommunicator( sub_ranks_tgt() );

    // Build the union communicator for the source and target. This is the
    // communicator over which we will operate the coupling.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_union;
    DataTransferKit::CommTools::unite( src_comm, tgt_comm, comm_union );

    // Set a boolean for source/target existence.
    bool src_exists = false;
    if ( !src_comm.is_null() ) src_exists = true;
    bool tgt_exists = false;
    if ( !tgt_comm.is_null() ) tgt_exists = true;
    comm_union->barrier();

    // Create a local to global process indexer for the source
    DataTransferKit::CommIndexer src_indexer( comm_union, src_comm );



    // ---------------//
    // SOURCE SETUP
    // ---------------//

    // Set required variables in the scope of the global communicator.
    Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Box,int> >
	src_geometry;
    Teuchos::RCP<
	DataTransferKit::FieldEvaluator<int,DataTransferKit::FieldContainer<double> > >
	src_evaluator;

    // If the source code exists on this process, build its data structures.
    if ( src_exists )
    {
	// Get source geometry.
	src_geometry = createSourceGeometry( src_comm );

	// Get the source field evaluator.
	src_evaluator = Teuchos::rcp( 
	    new SourceEvaluator( src_geometry->gids() ) );
    }
    comm_union->barrier();



    // ---------------//
    // TARGET SETUP
    // ---------------//

    // Set required variables in the scope of the global communicator.
    Teuchos::RCP<DataTransferKit::FieldManager<DataTransferKit::MeshContainer<int> > >
	tgt_coords;
    Teuchos::ArrayRCP<double> tgt_array;
    Teuchos::RCP<
	DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> > > 
	tgt_space;

    // If the target code exists on this process, build its data structures.
    if ( tgt_exists )
    {
	// Get the target mesh.
	Teuchos::RCP<DataTransferKit::MeshContainer<int> > tgt_mesh = 
	    createTargetMesh( tgt_comm );

	// Extract the mesh vertices as the target coordinates.
	tgt_coords = Teuchos::rcp( 
	    new DataTransferKit::FieldManager<DataTransferKit::MeshContainer<int> >(
		tgt_mesh, tgt_comm ) );

	// Get the target data space.
	tgt_array = Teuchos::ArrayRCP<double>( 12 );
	Teuchos::RCP<DataTransferKit::FieldContainer<double> > tgt_container = 
	    Teuchos::rcp( new DataTransferKit::FieldContainer<double>( tgt_array, 1 ) );
	tgt_space = Teuchos::rcp(
	    new DataTransferKit::FieldManager<DataTransferKit::FieldContainer<double> >(
		tgt_container, tgt_comm ) );
    }
    comm_union->barrier();
    


    // ---------------//
    // MAPPING SETUP
    // ---------------//

    // Create the mapping for the source-to-target transfer. The mapping will
    // occur over the union communicator in 3 dimensions. Keep track of missed
    // points as we expect to miss some.
    DataTransferKit::VolumeSourceMap<DataTransferKit::Box,
				     int,
				     DataTransferKit::MeshContainer<int> >
	src_to_tgt_map( comm_union, 3, true );

    // Setup the source-to-target map with the source geometry as the source
    // and the target coordinates as the target.
    src_to_tgt_map.setup( src_geometry, tgt_coords );



    // ---------------//
    // DATA TRANSFER
    // ---------------//
    
    // Apply the mapping to transfer the data.
    src_to_tgt_map.apply( src_evaluator, tgt_space );

    Teuchos::ArrayView<int> missed = src_to_tgt_map.getMissedTargetPoints();

    // Check the resulting data transfer.
    int fail_count = 0;
    if ( tgt_exists )
    {
	int target_rank = tgt_comm->getRank();

	// First, check for missed points. 
	// Proc 0 in the target decomposition should have 0.
	if ( target_rank == 0 )
	{
	    if ( missed.size() > 0 )
	    {
		++fail_count;
	    }
	}
	// Proc 1 in the target decomposition should have 4.
	else if ( target_rank == 1 )
	{
	    if ( missed.size() != 4 )
	    {
		++fail_count;
	    }
	}

	// Next, check the data. 
	if ( target_rank == 0 )
	{
	    // Proc 0 gets 1.0 or 2.0 or 3.0 and no misses.
	    for ( int i = 0; i < 12; ++i )
	    {
		if ( tgt_array[i] != 1.0 &&
		     tgt_array[i] != 2.0 &&
		     tgt_array[i] != 3.0 )
		{
		    ++fail_count;
		}
	    }
	}
	else if ( target_rank == 1 )
	{
	    for ( int i = 0; i < 12; ++i )
	    {
		// If we missed the point, should have got 0.
		if ( std::binary_search( missed.begin(), missed.end(), i ) )
		{
		    if ( tgt_array[i] != 0.0 )
		    {
			++fail_count;
		    }
		}
		// Otherwise we should have got 2.0 or 3.0.
		else if ( tgt_array[i] != 2.0 &&
			  tgt_array[i] != 3.0 )
		{
		    ++fail_count;
		}
	    }

	}
		     

	std::cout << std::endl;
	if ( fail_count == 0 )
	{
	    std::cout << "TEST PASSED: target proc " 
		      << tgt_comm->getRank() << std::endl;
	}
	else
	{
	    std::cout << "TEST FAILED " << fail_count << ": target proc " 
		      << tgt_comm->getRank() << std::endl;
	}
	std::cout << std::endl;
    }
    comm_union->barrier();

    return 0;
}

//---------------------------------------------------------------------------//
// end cxx_main.hpp
//---------------------------------------------------------------------------//

