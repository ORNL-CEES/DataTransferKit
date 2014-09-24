//---------------------------------------------------------------------------//
/*!
 * \file   DTK_DBC.cpp
 * \author Stuart Slattery
 * \brief  Assertions for error handling and Design-by-Contract.
 */
//---------------------------------------------------------------------------//

#include <sstream>

#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Assertion functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Build an assertion output from advanced constructor arguments.
 *
 * \param cond A string containing the assertion condition that failed.
 *
 * \param field A string containing the file name in which the assertion
 * failed. 
 *
 * \param line The line number at which the assertion failed.
 *
 * \return Assertion output.
 */
std::string Assertion::generate_output( 
    const std::string& cond, const std::string& file, const int line ) const
{
    std::ostringstream output;
    output << "DataTransferKit Assertion: " << cond << ", failed in " << file
	   << ", line " << line  << "." << std::endl;
    return output.str();
}

//---------------------------------------------------------------------------//
// Throw functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Throw a DataTransferKit::Assertion.
 *
 * \param cond A string containing the assertion condition that failed.
 *
 * \param field A string containing the file name in which the assertion
 * failed. 
 *
 * \param line The line number at which the assertion failed.
 */
void throwAssertion( const std::string& cond, const std::string& file,
		     const int line )
{
    throw Assertion( cond, file, line );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Throw a DataTransferKit::Assertion when an error code fails.
 *
 * \param cond A string containing the assertion condition that failed.
 *
 * \param field A string containing the file name in which the assertion
 * failed. 
 *
 * \param line The line number at which the assertion failed.
 *
 * \param error_code
 */
void errorCodeFailure( const std::string& cond, const std::string& file,
		       const int line, const int error_code )
{
    std::ostringstream output_msg;
    output_msg <<  "Error code : " << cond << ", failed in "
	      << file << ":" << line << std::endl
	      << "with error code:" << std::endl
	      << "\"" << error_code << "\"" << std::endl;
    throw Assertion( output_msg.str() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_DBC.cpp
//---------------------------------------------------------------------------//
