#include <iostream>

#include "Wave.hpp"
#include "Wave_Source.hpp"
#include "Wave_Target.hpp"
#include "Damper.hpp"
#include "Damper_Source.hpp"
#include "Damper_Target.hpp"

#include <DataTransferKit_DataSource.hpp>
#include <DataTransferKit_DataTarget.hpp>
#include <DataTransferKit_CopyOperator.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

//---------------------------------------------------------------------------//
// Main function driver for the coupled Wave/Damper problem.
int main(int argc, char* argv[])
{
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

    // Setup a Wave.
    Teuchos::RCP<Wave> wave =
	Teuchos::rcp( new Wave(comm, myMin, myMax, 10) );

    // Setup a Damper.
    Teuchos::RCP<Damper> damper =
	Teuchos::rcp( new Damper(comm, myMin, myMax, 10) ); 

    // Setup a Wave Data Source for the wave field.
    Teuchos::RCP<DataTransferKit::DataSource<double,int,double,1> > 
	wave_source = Teuchos::rcp( 
	    new DataTransferKit::Wave_DataSource<double,int,double,1>(wave) );

    // Setup a Damper Data Target for the wave field.
    Teuchos::RCP<DataTransferKit::DataTarget<double,int,double,1> > 
	damper_target = Teuchos::rcp( 
	    new DataTransferKit::Damper_DataTarget<double,int,double,1>(damper) );

    // Setup a copy operator for the wave field.
    DataTransferKit::CopyOperator<double,int,double,1> 
	wave_field_op( comm,
		       "WAVE_SOURCE_FIELD",
		       "WAVE_TARGET_FIELD",
		       wave_source,
		       damper_target );

    // Setup a Damper Data Source for the damper field.
    Teuchos::RCP<DataTransferKit::DataSource<double,int,double,1> > 
	damper_source = Teuchos::rcp( 
	    new DataTransferKit::Damper_DataSource<double,int,double,1>(damper) );

    // Setup a Wave Data Target for the damper field.
    Teuchos::RCP<DataTransferKit::DataTarget<double, int, double,1> > 
	wave_target = Teuchos::rcp( 
	    new DataTransferKit::Wave_DataTarget<double,int,double,1>(wave) );

    // Setup a copy operator for the damper field.
    DataTransferKit::CopyOperator<double,int,double,1> 
	damper_field_op( comm,
			 "DAMPER_SOURCE_FIELD",
			 "DAMPER_TARGET_FIELD",
			 damper_source,
			 wave_target );

    // Create the mapping for the wave field.
    wave_field_op.create_copy_mapping();

    // Create the mapping for the damper field.
    damper_field_op.create_copy_mapping();

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

	// Collect the norm values from the wave solve to check
	// convergence. 
	Teuchos::reduceAll<int>( *comm,
				 Teuchos::REDUCE_MAX, 
				 int(1), 
				 &local_norm, 
				 &global_norm );

	// Update the iteration count.
	++num_iter;

	// Barrier before proceeding.
	Teuchos::barrier<int>( *comm );
    }

    // Output results.
    if ( myRank == 0 )
    {
	std::cout << "Iterations to converge: " << num_iter << std::endl;
	std::cout << "L2 norm:                " << global_norm << std::endl;
    }

    return 0;
}
