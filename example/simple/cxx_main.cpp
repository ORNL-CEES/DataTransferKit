#include <algorithm>

#include <Wave.hpp>
#include <Damper.hpp>

#include <Coupler_Data_Source.hpp>
#include <Coupler_Data_Target.hpp>
#include <Coupler_Data_Field.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

//---------------------------------------------------------------------------//
// Get the current default communicator.
template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Main function driver for the coupled Wave/Damper problem.
int main(int argc, char* argv[])
{
    // Setup MPI.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    // Setup the parallel domain.
    double myMin;
    double myMax;

    // Setup a Wave.
    Teuchos::RCP<Wave> wave = Teuchos::rcp( new Wave(myMin, myMax, 10) );

    // Setup a Damper.
    Teuchos::RCP<Damper> damper = Teuchos::rcp( new Damper(myMin, myMax, 10) );

    // Setup a Wave Data Source for the wave field.

    // Setup a Damper Data Target for the wave field.

    // Setup a Data Field for the wave field.

    // Setup a Damper Data Source for the damper field.

    // Setup a Wave Data Target for the damper field.

    // Setup a Data Field for the damper field.

    // Iterate between the damper and wave until convergence.
    double norm = 1.0;
    int num_iter = 0;
    int max_iter = 100;
    while( norm < 1.0e-6 && num_iter < max_iter )
    {
	// Transfer the wave field.
	wave_field.transfer();

	// Damper solve.
	damper.solve();

	// Transfer the damper field.
	damper_field.transfer();

	// Wave solve.
	norm = wave.solve();

	// Update the iteration count.
	++num_iter;
    }

    return 0;
}
