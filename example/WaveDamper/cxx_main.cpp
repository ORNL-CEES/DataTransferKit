//---------------------------------------------------------------------------//
/*!
 * \file cxx_main.cpp
 * \author Stuart R. Slattery
 * \brief Multi-physics driver for WaveDamper example.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "Wave.hpp"
#include "Damper.hpp"

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

    // Create a Wave
    Teuchos::RCP<Wave> wave =
	Teuchos::rcp( new Wave(comm, myMin, myMax, 10) );

    // Get the Wave mesh.

    // Get the Wave field evaluator.

    // Get the Wave target coordinates.

    // Get the Wave target space.


    // ---------------//
    // DAMPER SETUP
    // ---------------//

    // Create a Damper.
    Teuchos::RCP<Damper> damper =
	Teuchos::rcp( new Damper(comm, myMin, myMax, 10) ); 

    // Get the Damper mesh.

    // Get the Damper field evaluator.

    // Get the Damper target coordinates.

    // Get the Damper target space.

    
    // ---------------//
    // MAPPING SETUP
    // ---------------//

    // Build the union communicator for the Wave and Damper.

    // Create the mapping for the wave-to-damper transfer.

    // Create the mapping for the damper-to-wave transfer.


    // ---------------//
    // SOLVE
    // ---------------//

    // Iterate between the damper and wave until convergence.
    double local_norm = 0.0;
    double global_norm = 1.0;
    int num_iter = 0;
    int max_iter = 100;
    while( global_norm > 1.0e-6 && num_iter < max_iter )
    {
	// Transfer the wave field.
	wave_field_op.copy();

	// Damper solve.
	damper->solve();

	// Transfer the damper field.
	damper_field_op.copy();

	// Wave solve.
	local_norm = wave->solve();

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
	std::cout << "L2 norm:                " << global_norm << std::endl;
    }

    return 0;
}
