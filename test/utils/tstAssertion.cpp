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

#include <DataTransferKit_config.hpp>
#include <DTK_Assertion.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
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
	testPrecondition( 0 );
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
	testPostcondition( 0 );
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
	testInvariant( 0 );
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
    rememberValue( test_value_1 = 0 );
    rememberValue( test_value_2 = 1 );
 
    try 
    {
	testInvariant( test_value_1 );
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
	testInvariant( test_value_2 );
	TEST_ASSERT( 1 );
    }
    catch( ... )
    {
	TEST_ASSERT( 0 );
    }
}

//---------------------------------------------------------------------------//
// end tstAssertion.cpp
//---------------------------------------------------------------------------//
