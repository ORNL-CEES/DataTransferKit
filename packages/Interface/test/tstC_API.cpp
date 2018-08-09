#include <Kokkos_Core.hpp>

#include <DTK_C_API.hpp>
#include <DTK_UserApplication.hpp>

#include "DTK_TestApplicationHelpers.hpp"

#include <string>

template <class UserApplication>
int test( const std::string &test_name, UserApplication &user_app )
{
    bool success = true;
    try
    {
        Teuchos::RCP<Teuchos::FancyOStream> fancy =
            Teuchos::fancyOStream( Teuchos::rcpFromRef( std::cout ) );
        Teuchos::FancyOStream &out = *fancy;

        if ( test_name == "test_node_list" )
            test_node_list( user_app, out, success );
        else if ( test_name == "test_bounding_volume_list" )
            test_bounding_volume_list( user_app, out, success );
        else if ( test_name == "test_polyhedron_list" )
            test_polyhedron_list( user_app, out, success );
        else if ( test_name == "test_multiple_topology_cell" )
            test_multiple_topology_cell( user_app, out, success );
        else if ( test_name == "test_boundary" )
            test_boundary( user_app, out, success );
        else if ( test_name == "test_adjacency_list" )
            test_adjacency_list( user_app, out, success );
        else if ( test_name == "test_single_topology_dof" )
            test_single_topology_dof( user_app, out, success );
        else if ( test_name == "test_multiple_topology_dof" )
            test_multiple_topology_dof( user_app, out, success );
        else if ( test_name == "test_field_push_pull" )
            test_field_push_pull( user_app, out, success );
        else if ( test_name == "test_field_eval" )
            test_field_eval( user_app, out, success );
        else if ( test_name == "test_missing_function" )
            test_missing_function( user_app, out, success );
        else if ( test_name == "test_too_many_functions" )
            test_too_many_functions( user_app, out, success );
        else
            throw std::runtime_error( "Unknown test name" );
    }
    catch ( ... )
    {
        return 2;
    }

    return success ? 0 : 1;
}

extern "C" {

int check_registry( const char *name, DTK_UserApplicationHandle handle )
{
    using namespace DataTransferKit;

    std::string test_name( name );

    DTK_Registry *dtk = reinterpret_cast<DTK_Registry *>( handle );
    auto registry = dtk->_registry;

    switch ( dtk->_space )
    {
    case DTK_SERIAL:
#ifdef KOKKOS_HAVE_SERIAL
    {
        UserApplication<double, Kokkos::HostSpace> user_app( registry );
        return_val = test( test_name, user_app );
    }
#else
        std::cout << "Serial execution space is disabled" << std::endl;
        return EXIT_FAILURE;
#endif
    case DTK_OPENMP:
#ifdef KOKKOS_HAVE_OPENMP
    {
        UserApplication<double, Kokkos::CudaUVMSpace> user_app( registry );
        return_val = test( test_name, user_app );
    }
#else
        std::cout << "OpenMP execution space is disabled" << std::endl;
        return EXIT_FAILURE;
#endif
    case DTK_CUDA:
#ifdef KOKKOS_HAVE_CUDA
    {
        UserApplication<double, Kokkos::Cuda> user_app( registry );
        return test( test_name, user_app );
    }
#else
        std::cout << "Cuda execution space is disabled" << std::endl;
        return EXIT_FAILURE;
#endif
    }

    return EXIT_SUCCESS;
}

} // extern "C"
