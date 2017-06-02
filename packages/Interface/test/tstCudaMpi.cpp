#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include <Kokkos_Core.hpp>

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// This is a 2 processors test.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CudaMpi, 1d_view, DeviceType )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Teuchos::RCP<const Teuchos::Comm<int>> comm_default =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm_default->getRank();
    unsigned int size = 10;

    if (comm_rank == 0)
    {
      int device = 0;
      cudaSetDevice(device);
      Kokkos::View<int*, DeviceType> view_rank("view_rank_0",size);
      Kokkos::parallel_for( "fill_rank_0",
                            Kokkos::RangePolicy<ExecutionSpace> (0, size),
                            KOKKOS_LAMBDA(int i) { view_rank[i] = 0; });
      Kokkos::fence();

	    auto view_rank_host = Kokkos::create_mirror_view(view_rank);
	  	Kokkos::deep_copy(view_rank_host, view_rank);
	  	for ( unsigned int i = 0; i < size; ++i )
	  		TEST_EQUALITY( view_rank_host( i ), 0 );

	  	int source = 1;
	  	int tag = 0;
	  	MPI_Status status;
    	MPI_Recv(&view_rank[0], size, MPI_INT, source, tag, MPI_COMM_WORLD, &status);

	  	Kokkos::deep_copy(view_rank_host, view_rank);
	  	for ( unsigned int i = 0; i < size; ++i )
	  		TEST_EQUALITY( view_rank_host( i ), 1 );
    }
    else if ( comm_rank == 1 )
    {
      int device = 1;
      cudaSetDevice(device);
      Kokkos::View<int*, DeviceType> view_rank("view_rank_1",size);
      Kokkos::parallel_for( "fill_rank_1",
                            Kokkos::RangePolicy<ExecutionSpace> (0, size),
                            KOKKOS_LAMBDA(int i) { view_rank[i] = 1; });
      Kokkos::fence();

	  	int destination = 0;
	  	int tag = 0;
	  	MPI_Send(&view_rank[0], size, MPI_INT, destination, tag, MPI_COMM_WORLD);
    }
}


// Include the test macros.
#include "DataTransferKitInterface_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CudaMpi, 1d_view, DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
