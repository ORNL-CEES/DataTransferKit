//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file   tstAssertion.cpp
 * \author Stuart Slattery
 * \brief  Assertion class unit tests.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <DTK_Assertion.hpp>

#include "Teuchos_UnitTestHarness.hpp"

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that a DataTransferKit::Assertion looks different than a
// std::runtime_error as it inherits from std::logic_error.
TEUCHOS_UNIT_TEST( Assertion, differentiation_test )
{
    try
    {
        throw std::runtime_error( "runtime error" );
    }
    catch ( const DataTransferKit::Assertion &assertion )
    {
        TEST_ASSERT( 0 );
    }
    catch ( ... )
    {
        TEST_ASSERT( 1 );
    }
}

//---------------------------------------------------------------------------//
// Check that a DataTransferKit::Assertion can be caught and the appropriate
// error message is written.
TEUCHOS_UNIT_TEST( Assertion, message_test )
{
    std::string message;

    try
    {
        throw DataTransferKit::Assertion( "cond", "file", 12 );
    }
    catch ( const DataTransferKit::Assertion &assertion )
    {
        message = std::string( assertion.what() );
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }

    const std::string true_message(
        "DataTransferKit Assertion: cond, failed in file, line 12.\n" );
    TEST_ASSERT( 0 == message.compare( true_message ) );
}

//---------------------------------------------------------------------------//
// Check that we can throw a nemesis assertion with throwAssertion.
TEUCHOS_UNIT_TEST( Assertion, throw_test )
{
    try
    {
        const std::string message( "message" );
        const std::string file( "file" );
        const int line( 12 );
        DataTransferKit::throwAssertion( message, file, line );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::Assertion &assertion )
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
TEUCHOS_UNIT_TEST( Assertion, precondition_test )
{
    try
    {
        testPrecondition( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::Assertion &assertion )
    {
#if HAVE_DTK_DBC
        std::string message( assertion.what() );
        std::string true_message( "DataTransferKit Assertion: 0, failed in" );
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
TEUCHOS_UNIT_TEST( Assertion, postcondition_test )
{
    try
    {
        testPostcondition( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::Assertion &assertion )
    {
#if HAVE_DTK_DBC
        std::string message( assertion.what() );
        std::string true_message( "DataTransferKit Assertion: 0, failed in" );
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
TEUCHOS_UNIT_TEST( Assertion, invariant_test )
{
    try
    {
        testInvariant( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::Assertion &assertion )
    {
#if HAVE_DTK_DBC
        std::string message( assertion.what() );
        std::string true_message( "DataTransferKit Assertion: 0, failed in" );
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
TEUCHOS_UNIT_TEST( Assertion, remember_test )
{
    remember( int test_value_1 = 0 );
    remember( int test_value_2 = 1 );

    try
    {
        testInvariant( test_value_1 );
    }
    catch ( const DataTransferKit::Assertion &assertion )
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
        testInvariant( test_value_2 );
        TEST_ASSERT( 1 );
    }
    catch ( ... )
    {
        TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Test the assertion check for DBC.
TEUCHOS_UNIT_TEST( Assertion, assertion_test )
{
    try
    {
        insist( 0 );
        throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch ( const DataTransferKit::Assertion &assertion )
    {
        std::string message( assertion.what() );
        std::string true_message( "DataTransferKit Assertion: 0, failed in" );
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
// end tstAssertion.cpp
//---------------------------------------------------------------------------//
