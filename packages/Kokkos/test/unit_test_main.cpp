#include <Kokkos_Core.hpp>

#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_GlobalMPISession.hpp>


int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  Kokkos::initialize(argc, argv);
  int return_val = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  Kokkos::finalize();
  return return_val;
}
