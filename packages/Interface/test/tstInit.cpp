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

    auto check = [&status]( bool cond ) { status = ( status && cond ); };

    if ( t == 1 )
    {
        check( !DataTransferKit::isInitialized() &&
               !Kokkos::DefaultExecutionSpace::is_initialized() );

        DataTransferKit::initialize( &argc, &argv );
        check( DataTransferKit::isInitialized() &&
               Kokkos::DefaultExecutionSpace::is_initialized() );

        DataTransferKit::finalize();
        check( !DataTransferKit::isInitialized() &&
               !Kokkos::DefaultExecutionSpace::is_initialized() );
    }
    else if ( t == 2 )
    {
        Kokkos::initialize( argc, argv );

        check( !DataTransferKit::isInitialized() &&
               Kokkos::DefaultExecutionSpace::is_initialized() );

        DataTransferKit::initialize( &argc, &argv );
        check( DataTransferKit::isInitialized() &&
               Kokkos::DefaultExecutionSpace::is_initialized() );

        DataTransferKit::finalize();
        check( !DataTransferKit::isInitialized() &&
               Kokkos::DefaultExecutionSpace::is_initialized() );

        Kokkos::finalize();
        check( !Kokkos::DefaultExecutionSpace::is_initialized() );
    }
    else if ( t == 3 )
    {
        check( !DataTransferKit::isInitialized() &&
               !Kokkos::DefaultExecutionSpace::is_initialized() );

        DataTransferKit::initialize();
        check( DataTransferKit::isInitialized() &&
               Kokkos::DefaultExecutionSpace::is_initialized() );

        DataTransferKit::finalize();
        check( !DataTransferKit::isInitialized() &&
               !Kokkos::DefaultExecutionSpace::is_initialized() );
    }
    else
    {
        status = false;
    }

    std::cout << "End Result: " << ( status ? "TEST PASSED" : "TEST FAILED" )
              << std::endl;
    return 0;
}
