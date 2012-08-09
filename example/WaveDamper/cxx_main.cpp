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

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>

//---------------------------------------------------------------------------//
// Main function driver for the coupled Wave/Damper problem.
int main(int argc, char* argv[])
{
    // ---------------//
    // PARALLEL SETUP
    // ---------------//

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Set up the parallel domain.
    double global_min = 0.0;
    double global_max = 5.0;
    int myRank = comm->getRank();
    int mySize = comm->getSize();
    double local_size = (global_max - global_min) / mySize;
    double myMin = myRank*local_size + global_min;
    double myMax = (myRank+1)*local_size + global_min;


    // ---------------//
    // WAVE SETUP
    // ---------------//

    // Create a Wave.
    Teuchos::RCP<Wave> wave =
	Teuchos::rcp( new Wave(comm, myMin, myMax, 10) );

    // Get the Wave mesh.
    Teuchos::RCP<DataTransferKit::MeshManager<WaveAdapter::MeshType> > 
	wave_mesh = WaveAdapter::getMesh( wave );

    // Get the Wave field evaluator.
    WaveAdapter::RCP_Evaluator wave_evaluator =
	WaveAdapter::getFieldEvaluator( wave );

    // Get the Wave target coordinates.
    Teuchos::RCP<DataTransferKit::FieldManager<WaveAdapter::MeshType> >
	wave_target_coords = WaveAdapter::getTargetCoords( wave );

    // Get the Wave target space.
    Teuchos::RCP<DataTransferKit::FieldManager<WaveAdapter::FieldType> >
	wave_target_space = WaveAdapter::getTargetSpace( wave );


    // ---------------//
    // DAMPER SETUP
    // ---------------//

    // Create a Damper.
    Teuchos::RCP<Damper> damper =
	Teuchos::rcp( new Damper(comm, myMin, myMax, 10) ); 

    // Get the Damper mesh.
    Teuchos::RCP<DataTransferKit::MeshManager<DamperAdapter::MeshType> >
	damper_mesh = DamperAdapter::getMesh( damper );

    // Get the Damper field evaluator.
    DamperAdapter::RCP_Evaluator damper_evaluator =
	DamperAdapter::getFieldEvaluator( damper );

    // Get the Damper target coordinates.
    Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::MeshType> >
	damper_target_coords = DamperAdapter::getTargetCoords( damper );

    // Get the Damper target space.
    Teuchos::RCP<DataTransferKit::FieldManager<DamperAdapter::FieldType> >
	damper_target_space = DamperAdapter::getTargetSpace( damper );

    
    // ---------------//
    // MAPPING SETUP
    // ---------------//

    // Build the union communicator for the Wave and Damper.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_union;
    DataTransferKit::CommTools::unite( wave->get_comm(), damper->get_comm(),
				       comm_union );

    // Create the mapping for the wave-to-damper transfer.
    DataTransferKit::SharedDomainMap<WaveAdapter::MeshType,DamperAdapter::MeshType>
	wave_to_damper_map( comm_union );
    wave_to_damper_map.setup( wave_mesh, damper_target_coords );

    // Create the mapping for the damper-to-wave transfer.
    DataTransferKit::SharedDomainMap<DamperAdapter::MeshType,WaveAdapter::MeshType>
	damper_to_wave_map( comm_union );
    damper_to_wave_map.setup( damper_mesh, wave_target_coords );


    // ---------------//
    // SOLVE
    // ---------------//

    // Iterate between the damper and wave until convergence.
    double norm = 0.0;
    int num_iter = 0;
    int max_iter = 100;
    while( norm > 1.0e-6 && num_iter < max_iter )
    {
	// Transfer the wave field.
	wave_to_damper_map.apply( wave_evaluator, damper_target_space );

	// Damper solve.
	damper->solve();

	// Transfer the damper field.
	damper_to_wave_map.apply( damper_evaluator, wave_target_space );

	// Wave solve.
	wave->solve();

	// Get the norm of the wave field to check convergence.
	

	// Update the iteration count.
	++num_iter;

	// Barrier before proceeding.
	comm->barrier();
    }

    // Output results.
    if ( myRank == 0 )
    {
	std::cout << "Iterations to converge: " << num_iter << std::endl;
	std::cout << "L2 norm:                " << norm << std::endl;
    }

    return 0;
}

//---------------------------------------------------------------------------//
// end cxx_main.hpp
//---------------------------------------------------------------------------//

