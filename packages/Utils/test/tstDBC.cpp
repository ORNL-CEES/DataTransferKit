/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file   tstDBC.cpp
 * \author Stuart Slattery
 * \brief  Assertion class unit tests.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <DTK_DBC.hpp>

#include "Teuchos_UnitTestHarness.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
int error_code_function( int &value )
{
    int code = value;
    ++value;
    return code;
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that a DataTransferKit::DataTransferKitException looks different than a
// std::runtime_error as it inherits from std::logic_error.
TEUCHOS_UNIT_TEST( DataTransferKitException, differentiation_test )
{
    try
    {
        throw std::runtime_error( "runtime error" );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
        TEST_ASSERT( 0 );
    }
    catch ( ... )
    {
        TEST_ASSERT( 1 );
    }
}

//---------------------------------------------------------------------------//
// Check that a DataTransferKit::DataTransferKitException can be caught and the
// appropriate
// error message is written.
TEUCHOS_UNIT_TEST( DataTransferKitException, message_test )
{
    std::string message;

    try
    {
        throw DataTransferKit::DataTransferKitException( "cond", "file", 12 );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
        message = std::string( assertion.what() );
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }

    const std::string true_message( "DataTransferKit DataTransferKitException: "
                                    "cond, failed in file, line 12.\n" );
    TEST_ASSERT( 0 == message.compare( true_message ) );
}

//---------------------------------------------------------------------------//
// Check that we can throw a nemesis assertion with
// throwDataTransferKitException.
TEUCHOS_UNIT_TEST( DataTransferKitException, throw_test )
{
    try
    {
        const std::string message( "message" );
        const std::string file( "file" );
        const int line( 12 );
        DataTransferKit::throwDataTransferKitException( message, file, line );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
        TEST_ASSERT( 1 );
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Test the precondition check for DBC.
TEUCHOS_UNIT_TEST( DataTransferKitException, precondition_test )
{
    try
    {
        DTK_REQUIRE( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
#if HAVE_DTK_DBC
        std::string message( assertion.what() );
        std::string true_message(
            "DataTransferKit DataTransferKitException: 0, failed in" );
        std::string::size_type idx = message.find( true_message );
        if ( idx == std::string::npos )
        {
            TEST_ASSERT( 0 );
        }
#else
        TEST_ASSERT( 0 );
#endif
    }
    catch ( ... )
    {
#if HAVE_DTK_DBC
        TEST_ASSERT( 0 );
#endif
    }
}

//---------------------------------------------------------------------------//
// Test the postcondition check for DBC.
TEUCHOS_UNIT_TEST( DataTransferKitException, postcondition_test )
{
    try
    {
        DTK_ENSURE( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
#if HAVE_DTK_DBC
        std::string message( assertion.what() );
        std::string true_message(
            "DataTransferKit DataTransferKitException: 0, failed in" );
        std::string::size_type idx = message.find( true_message );
        if ( idx == std::string::npos )
        {
            TEST_ASSERT( 0 );
        }
#else
        TEST_ASSERT( 0 );
#endif
    }
    catch ( ... )
    {
#if HAVE_DTK_DBC
        TEST_ASSERT( 0 );
#endif
    }
}

//---------------------------------------------------------------------------//
// Test the invariant check for DBC.
TEUCHOS_UNIT_TEST( DataTransferKitException, invariant_test )
{
    try
    {
        DTK_CHECK( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
#if HAVE_DTK_DBC
        std::string message( assertion.what() );
        std::string true_message(
            "DataTransferKit DataTransferKitException: 0, failed in" );
        std::string::size_type idx = message.find( true_message );
        if ( idx == std::string::npos )
        {
            TEST_ASSERT( 0 );
        }
#else
        TEST_ASSERT( 0 );
#endif
    }
    catch ( ... )
    {
#if HAVE_DTK_DBC
        TEST_ASSERT( 0 );
#endif
    }
}

//---------------------------------------------------------------------------//
// Test that we can remember a value and check it with DBC.
TEUCHOS_UNIT_TEST( DataTransferKitException, remember_test )
{
    DTK_REMEMBER( int test_value_1 = 0 );
    DTK_REMEMBER( int test_value_2 = 1 );

    try
    {
        DTK_CHECK( test_value_1 );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
#if HAVE_DTK_DBC
        TEST_ASSERT( 1 );
#else
        TEST_ASSERT( 0 );
#endif
    }
    catch ( ... )
    {
#if HAVE_DTK_DBC
        TEST_ASSERT( 0 );
#endif
    }

    try
    {
        DTK_CHECK( test_value_2 );
        TEST_ASSERT( 1 );
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Test the assertion check for DBC.
TEUCHOS_UNIT_TEST( DataTransferKitException, assertion_test )
{
    try
    {
        DTK_INSIST( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
        std::string message( assertion.what() );
        std::string true_message(
            "DataTransferKit DataTransferKitException: 0, failed in" );
        std::string::size_type idx = message.find( true_message );
        if ( idx == std::string::npos )
        {
            TEST_ASSERT( 0 );
        }
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Test the error code check.
TEUCHOS_UNIT_TEST( DataTransferKitException, errorcode_test_1 )
{
    int value = 1;
    try
    {
        DTK_CHECK_ERROR_CODE( error_code_function( value ) );
        TEST_EQUALITY( value, 2 );
    }
    catch ( const DataTransferKit::DataTransferKitException &assertion )
    {
#if HAVE_DTK_DBC
        TEST_EQUALITY( value, 2 );
#else
        TEST_ASSERT( 0 );
#endif
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Test the error code check.
TEUCHOS_UNIT_TEST( DataTransferKitException, errorcode_test_2 )
{
    int value = 0;
    try
    {
        DTK_CHECK_ERROR_CODE( error_code_function( value ) );
        TEST_EQUALITY( value, 1 );
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Check that we can catch a
// DataTransferKit::DataTransferKitNotImplementedException and that the
// appropriate error message is written.
TEUCHOS_UNIT_TEST( DataTransferKitNotImplementedException, message_error )
{
    std::string message;
    try
    {
        throw DataTransferKit::DataTransferKitNotImplementedException();
    }
    catch ( const DataTransferKit::DataTransferKitNotImplementedException
                &assertion )
    {
        message = std::string( assertion.what() );
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }

    const std::string true_message( "Not implemented" );
    TEST_ASSERT( 0 == message.compare( true_message ) );
}

//---------------------------------------------------------------------------//
// end tstDBC.cpp
//---------------------------------------------------------------------------//
