//---------------------------------------------------------------------------//
/*!
 * \file   tstAssertion.cpp
 * \author Stuart Slattery
 * \brief  Assertion class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_DBC.hpp>

#include "Teuchos_UnitTestHarness.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
int error_code_function( int& value )
{
    int code = value;
    ++value;
    return code;
}

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
    catch( const DataTransferKit::Assertion& assertion )
    {
	TEST_ASSERT( 0 );
    }
    catch( ... )
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
    catch( const DataTransferKit::Assertion& assertion )
    {
	message = std::string( assertion.what() );
    }
    catch( ... )
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
    catch( const DataTransferKit::Assertion& assertion )
    {
	TEST_ASSERT( 1 );	
    }
    catch( ... )
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
	DTK_REQUIRE( 0 );
	throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch( const DataTransferKit::Assertion& assertion )
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
    catch( ... )
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
	DTK_ENSURE( 0 );
	throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch( const DataTransferKit::Assertion& assertion )
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
    catch( ... )
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
	DTK_CHECK( 0 );
	throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch( const DataTransferKit::Assertion& assertion )
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
    catch( ... )
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
    DTK_REMEMBER( int test_value_1 = 0 );
    DTK_REMEMBER( int test_value_2 = 1 );
 
    try 
    {
	DTK_CHECK( test_value_1 );
    }
    catch( const DataTransferKit::Assertion& assertion )
    {
#if HAVE_DTK_DBC
	TEST_ASSERT( 1 );
#else
	TEST_ASSERT( 0 );
#endif
    }
    catch( ... )
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
    catch( ... )
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
	DTK_INSIST( 0 );
	throw std::runtime_error( "this shouldn't be thrown" );
    }
    catch( const DataTransferKit::Assertion& assertion )
    {
	std::string message( assertion.what() );
	std::string true_message( "DataTransferKit Assertion: 0, failed in" );
	std::string::size_type idx = message.find( true_message );
	if ( idx == std::string::npos )
	{
	    TEST_ASSERT( 0 );
	}
    }
    catch( ... )
    {
	TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Test the error code check.
TEUCHOS_UNIT_TEST( Assertion, errorcode_test_1 )
{
    int value = 1;
    try 
    {
	DTK_CHECK_ERROR_CODE( error_code_function(value) );
	TEST_EQUALITY( value, 2 );
    }
    catch( const DataTransferKit::Assertion& assertion )
    {
#if HAVE_DTK_DBC
	TEST_EQUALITY( value, 2 );
#else
	TEST_ASSERT( 0 );
#endif
    }
    catch( ... )
    {
	TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// Test the error code check.
TEUCHOS_UNIT_TEST( Assertion, errorcode_test_2 )
{
    int value = 0;
    try 
    {
	DTK_CHECK_ERROR_CODE( error_code_function(value) );
	TEST_EQUALITY( value, 1 );
    }
    catch( ... )
    {
	TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// end tstAssertion.cpp
//---------------------------------------------------------------------------//
