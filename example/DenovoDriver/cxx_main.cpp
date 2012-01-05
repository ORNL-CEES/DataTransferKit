#include <iostream>

#include "DenovoSource.hpp"
#include "DenovoTarget.hpp"

#include <Coupler_DataSource.hpp>
#include <Coupler_DataTarget.hpp>
#include <Coupler_DataField.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

//---------------------------------------------------------------------------//
// Main function driver for the DenovoDriver example.
int main(int argc, char* argv[])
{
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();



    return 0;
}
