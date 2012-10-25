//---------------------------------------------------------------------------//
/*!
 * \file cxx_main.cpp
 * \author Stuart R. Slattery
 * \brief Multi-physics driver for WaveDamper example.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "Wave.hpp"
#include "WaveAdapter.hpp"
#include "Damper.hpp"
#include "DamperAdapter.hpp"

#include <DTK_SharedDomainMap.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_FieldTools.hpp>
#include <DTK_CommTools.hpp>
#include <DTK_CommIndexer.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_Ptr.hpp>

//---------------------------------------------------------------------------//
// Main function driver for the coupled Wave/Damper problem.
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

    // Split the main communicator into 2 separate groups, one from the wave
    // and one for the damper.
    Teuchos::Array<int> sub_ranks_wave, sub_ranks_damper;

    // Multiple processor case.
    if ( num_procs > 1 )
    {
	for ( int n = 0; n < num_procs; ++n )
	{
	    if ( n % 2 == 0 )
	    {
		sub_ranks_wave.push_back(n);
	    }
	    else
	    {
		sub_ranks_damper.push_back(n);
	    }
	}
    }
    // Single processor case.
    else
    {
	sub_ranks_wave.push_back(0);
	sub_ranks_damper.push_back(0);
    }

    // Generate the wave and damper communicators from the sub ranks.
    Teuchos::RCP<const Teuchos::Comm<int> > wave_comm = 
	comm_default->createSubcommunicator( sub_ranks_wave() );

    Teuchos::RCP<const Teuchos::Comm<int> > damper_comm = 
	comm_default->createSubcommunicator( sub_ranks_damper() );

    // Build the union communicator for the Wave and Damper. This is the
    // communicator over which we will operate the coupling.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_union;
    DataTransferKit::CommTools::unite( wave_comm, damper_comm, comm_union );

    // Set a boolean for Wave/Damper existence.
    bool wave_exists = false;
    if ( !wave_comm.is_null() ) wave_exists = true;
    bool damper_exists = false;
    if ( !damper_comm.is_null() ) damper_exists = true;
    comm_union->barrier();

    // Create a local to global process indexer for the Wave.
    DataTransferKit::CommIndexer wave_indexer( comm_union, wave_comm );

    // Global domain.
    double global_min = 0.0;
    double global_max = 5.0;


    // ---------------//
    // WAVE SETUP
    // ---------------//

    // Set required variables in the scope of the global communicator.
    Teuchos::RCP<Wave> wave;
    Teuchos::RCP<DataTransferKit::MeshManager<WaveAdapter::MeshType> > 
	wave_mesh;
    WaveAdapter::RCP_Evaluator wave_evaluator;
    Teuchos::RCP<DataTransferKit::FieldManager<WaveAdapter::MeshType> > 
	wave_target_coords;
    Teuchos::RCP<DataTransferKit::FieldManager<WaveAdapter::FieldType> > 
	wave_target_space;

    // If the wave code exists on this process, build its data structures.
    if ( wave_exists )
    {
	// Set up the wave parallel domain.
	int my_wave_rank = wave_comm->getRank();
	int my_wave_size = wave_comm->getSize();
	double wave_local_size = (global_max - global_min) / my_wave_size;
	double my_wave_min = my_wave_rank*wave_local_size + global_min;
	double my_wave_max = (my_wave_rank+1)*wave_local_size + global_min;
	int wave_num_local_elements = 11;

	// Create a Wave.
	wave = Teuchos::rcp( new Wave( wave_comm, my_wave_min, 
				       my_wave_max, wave_num_local_elements) );

	// Get the Wave mesh.
	wave_mesh = WaveAdapter::getMesh( wave );

	// Get the Wave field evaluator.
	wave_evaluator = WaveAdapter::getFieldEvaluator( wave );

	// Get the Wave target coordinates.
	wave_target_coords = WaveAdapter::getTargetCoords( wave );

	// Get the Wave target space.
	wave_target_space = WaveAdapter::getTargetSpace( wave );
    }
    comm_union->barrier();

    // ---------------//
    // DAMPER SETUP
    // ---------------//

    // Set required variables in the scope of the global communicator.
    Teuchos::RCP<Damper> damper;
    Teuchos::RCP<DataTransferKit::MeshManager<DamperAdapter::MeshType> > 
	damper_mesh;
    DamperAdapter::RCP_Evaluator damper_evaluator;
    Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::MeshType> > 
	damper_target_coords;
    Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::FieldType> > 
	damper_target_space;

    // If the damper code exists on this process, build its data structures.
    if ( damper_exists )
    {
	// Set up the damper parallel domain.
	int my_damper_rank = damper_comm->getRank();
	int my_damper_size = damper_comm->getSize();
	double damper_local_size = (global_max - global_min) / my_damper_size;
	double my_damper_min = my_damper_rank*damper_local_size + global_min;
	double my_damper_max = (my_damper_rank+1)*damper_local_size + global_min;
	int damper_num_local_elements = 11;
	
	// Create a Damper.
	damper = Teuchos::rcp( 
	    new Damper( damper_comm, my_damper_min, 
			my_damper_max, damper_num_local_elements) );

	// Get the Damper mesh.
	damper_mesh = DamperAdapter::getMesh( damper );

	// Get the Damper field evaluator.
	damper_evaluator = DamperAdapter::getFieldEvaluator( damper );

	// Get the Damper target coordinates.
	damper_target_coords = DamperAdapter::getTargetCoords( damper );

	// Get the Damper target space.
	damper_target_space = DamperAdapter::getTargetSpace( damper );
    }
    comm_union->barrier();
    
    // ---------------//
    // MAPPING SETUP
    // ---------------//

    // Create the mapping for the wave-to-damper transfer. The mapping will
    // occur over the union communicator in 1 dimensions.
    DataTransferKit::SharedDomainMap<WaveAdapter::MeshType,DamperAdapter::MeshType>
	wave_to_damper_map( comm_union, 1 );

    // Setup the wave-to-damper map with the wave mesh as the source and the
    // damper coordinates as the target.
    wave_to_damper_map.setup( wave_mesh, damper_target_coords );

    // Create the mapping for the damper-to-wave transfer. The mapping will
    // occur over the union communicator in 1 dimensions.
    DataTransferKit::SharedDomainMap<DamperAdapter::MeshType,WaveAdapter::MeshType>
	damper_to_wave_map( comm_union, 1 );

    // Setup the damper-to-wave map with the damper mesh as the source and the
    // wave coordinates as the target.
    damper_to_wave_map.setup( damper_mesh, wave_target_coords );


    // ---------------//
    // SOLVE
    // ---------------//

    // Setup a FieldContainer for the actual data in the Wave code so we can
    // compute its norm for the convergence check.
    Teuchos::Array<double> wave_norm( 1, 1.0 );
    Teuchos::RCP<std::vector<double> > wave_data = 
	Teuchos::rcp( new std::vector<double>(0,0) );
    if ( wave_exists ) wave_data = wave->get_data();
    comm_union->barrier();
    Teuchos::ArrayRCP<double> wave_arcp( &(*wave_data)[0], 0,
					 wave_data->size(), false );
    DataTransferKit::FieldContainer<double> wave_container( wave_arcp, 1 );

    // Iterate between the damper and wave until convergence.
    int num_iter = 0;
    int max_iter = 50;
    double norm = 1.0;
    double tolerance = 1.0e-6;
    while( norm > tolerance && num_iter < max_iter )
    {
	// Transfer the wave field.
	wave_to_damper_map.apply( wave_evaluator, damper_target_space );

	// Damper solve.
	if ( damper_exists ) 
	{
	    damper->solve();
	}
	comm_union->barrier();

	// Transfer the damper field.
	damper_to_wave_map.apply( damper_evaluator, wave_target_space );

	// Wave solve.
	if ( wave_exists ) 
	{
	    wave->solve();

	    DataTransferKit::FieldTools<
		DataTransferKit::FieldContainer<double> >::norm2(
		    wave_container, wave->get_comm(), wave_norm );

	    norm = wave_norm[0];
	}
	comm_union->barrier();
       
	// Broadcast the norm results from from Wave proc 0 to all processes
	// in the union communicator.
	Teuchos::broadcast<int,double>( *comm_union, wave_indexer.l2g(0),
					Teuchos::Ptr<double>(&norm) );

	// Update the iteration count.
	++num_iter;

	// Barrier before proceeding to the next iteration.
	comm_union->barrier();
    }

    // Output results from proc 0.
    if ( comm_union->getRank() == 0 )
    {
	std::cout << "Iterations to converge: " << num_iter << std::endl;
	std::cout << "L2 norm:                " << norm << std::endl;
    }
    comm_union->barrier();

    return 0;
}

//---------------------------------------------------------------------------//
// end cxx_main.hpp
//---------------------------------------------------------------------------//

