/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file
 * \brief  Kokkos initialization tests
 */

#include "DTK_Core.hpp"

#include <getopt.h>

int main( int argc, char *argv[] )
{
    bool status = true;

    // Parse command line parameters
    int opt;
    int t = 0;
    while ( ( opt = getopt( argc, argv, "t:h" ) ) != -1 )
    {
        switch ( opt )
        {
        case 't':
            t = atoi( optarg );
            break;

        case 'h':
            printf( "Usage: %s [-t <test_number>] [-h]\n", argv[0] );
            return EXIT_FAILURE;
        }
    }

    using default_space = Kokkos::DefaultExecutionSpace;
#ifdef KOKKOS_HAVE_SERIAL
    // Serial is always initialized when enabled. If it's the only execution
    // space, it's going to be the default, and thus going to be initialized
    // even without calling Kokkos::initialize(). Therefore, we need to skip
    // checking its initialization in this case.
    const bool kokkos_always_initialized =
        ( typeid( default_space ) == typeid( Kokkos::Serial ) );
#else
    const bool kokkos_always_initialized = false;
#endif

    auto check = [&status]( bool cond ) { status = ( status && cond ); };

    if ( t == 1 )
    {
        check(
            !DataTransferKit::isInitialized() &&
            ( kokkos_always_initialized || !default_space::is_initialized() ) );

        DataTransferKit::initialize( &argc, &argv );
        check( DataTransferKit::isInitialized() &&
               default_space::is_initialized() );

        DataTransferKit::finalize();
        check(
            !DataTransferKit::isInitialized() &&
            ( kokkos_always_initialized || !default_space::is_initialized() ) );
    }
    else if ( t == 2 )
    {
        Kokkos::initialize( argc, argv );
        check( !DataTransferKit::isInitialized() &&
               default_space::is_initialized() );

        DataTransferKit::initialize( &argc, &argv );
        check( DataTransferKit::isInitialized() &&
               default_space::is_initialized() );

        DataTransferKit::finalize();
        check( !DataTransferKit::isInitialized() &&
               default_space::is_initialized() );

        Kokkos::finalize();
        check( kokkos_always_initialized || !default_space::is_initialized() );
    }
    else if ( t == 3 )
    {
        check(
            !DataTransferKit::isInitialized() &&
            ( kokkos_always_initialized || !default_space::is_initialized() ) );

        DataTransferKit::initialize();
        check( DataTransferKit::isInitialized() &&
               default_space::is_initialized() );

        DataTransferKit::finalize();
        check(
            !DataTransferKit::isInitialized() &&
            ( kokkos_always_initialized || !default_space::is_initialized() ) );
    }
    else
    {
        status = false;
    }

    std::cout << "End Result: " << ( status ? "TEST PASSED" : "TEST FAILED" )
              << std::endl;
    return 0;
}
