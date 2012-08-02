//---------------------------------------------------------------------------//
/*!
 * \file   DataTransferKit_Exception.cpp
 * \author Stuart Slattery
 * \brief  Exceptions for error handling.
 */
//---------------------------------------------------------------------------//

#include <sstream>

#include "DTK_Exception.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Exception functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Build an exception output from advanced constructor arguments.
 */
std::string Exception::generate_output( 
    const std::string& cond, const std::string& file, const int line ) const
{
    std::ostringstream output;
    output << "DTK Exception: " << cond << ", failed in " << file 
	   << ", line " << line  << "." << std::endl;
    return output.str();
}

//---------------------------------------------------------------------------//
// Throw functions.
//---------------------------------------------------------------------------//
// Throw a DataTransferKit::Exception.
void throwException( const std::string& cond, const std::string& file,
		     const int line )
{
    throw Exception( cond, file, line );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Exception.cpp
//---------------------------------------------------------------------------//
