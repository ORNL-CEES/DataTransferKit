//---------------------------------------------------------------------------//
/*!
 * \file   DataTransferKit_Exception.hpp
 * \author Stuart Slattery
 * \brief  Exceptions for error handling.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EXCEPTION_HPP
#define DTK_EXCEPTION_HPP

#include <stdexcept>
#include <string>

namespace DataTransferKit
{

class PreconditionException : public std::runtime_error
{
  public:
    PreconditionException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

class PostconditionException : public std::runtime_error
{
  public:
    PostconditionException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

class InvariantException : public std::runtime_error
{
  public:
    InvariantException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

} // end namespace DataTransferKit

#endif // end DTK_EXCEPTION_HPP
