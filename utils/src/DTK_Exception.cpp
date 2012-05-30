//---------------------------------------------------------------------------//
/*!
 * \file   DataTransferKit_Exception.cpp
 * \author Stuart Slattery
 * \brief  Exceptions for error handling.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "DTK_Exception.hpp"

#include <Teuchos_TestForException.hpp>

namespace DataTransferKit
{
/*!
 * \brief Test for a precondition exception.
 */
void testPrecondition( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				PreconditionException,
				msg << std::endl );
}

/*!
 * \brief Test for a postcondition exception.
 */
void testPostcondition( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				PostconditionException,
				msg << std::endl );
}

/*!
 * \brief Test for a Invariant exception.
 */
void testInvariant( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				InvariantException,
				msg << std::endl );
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Exception.cpp
//---------------------------------------------------------------------------//
