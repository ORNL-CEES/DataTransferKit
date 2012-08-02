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
//---------------------------------------------------------------------------//
/*!
 * \brief Base class for DTK exceptions. This structure is heavily based on
 * that developed by Tom Evans.
 */
//---------------------------------------------------------------------------//
class Exception : public std::runtime_error
{
  public:

    /*! 
     * \brief Default constructor.
     */
    Exception( const std::string& msg )
	: std::runtime_error( msg )
    { /* ... */ }

    /*! 
     * \brief Advanced constructor.
     */
    Exception( const std::string& cond, const std::string& file, 
		  const int line )
	: std::runtime_error( generate_output( cond, file, line ) )
    { /* ... */ }

    //! Destructor.
    virtual ~Exception() throw()
    { /* ... */ }

  private:

    // Build an exception output from advanced constructor arguments.
    std::string generate_output( const std::string& cond, 
				 const std::string& file, 
				 const int line ) const;
};

//---------------------------------------------------------------------------//
// Throw functions.
//---------------------------------------------------------------------------//
// Throw a DataTransferKit::Exception.
void throwException( const std::string& cond, const std::string& file,
		     const int line );

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Design-by-Contract macros.
//---------------------------------------------------------------------------//
/*!
 * \page DataTransferKit Design-by-Contract.
 *
 * Design-by-Contract functionality is provided to verify function
 * preconditions, postconditions, and invariants. These checks are separated
 * from the debug build by setting the following in a CMake configure:
 *
 * -D DataTransferKit_ENABLE_DBC:BOOL=ON
 *
 * Although they will require computational overhead, these checks provide an
 * initial mechanism for veryifing library input arguments. Note that the
 * bounds-checking functionality used within the DataTransferKit is only
 * provided by a debug build.
 */

#if HAVE_DTK_DBC

#define testPrecondition(c) \
    if (!(c)) DataTransferKit::throwException( #c, __FILE__, __LINE__ )
#define testPostcondition(c) \
    if (!(c)) DataTransferKit::throwException( #c, __FILE__, __LINE__ )
#define testInvariant(c) \
    if (!(c)) DataTransferKit::throwException( #c, __FILE__, __LINE__ )

#else

#define testPrecondition(c)
#define testPostcondition(c)
#define testInvariant(c)

#endif

//---------------------------------------------------------------------------//

#endif // end DTK_EXCEPTION_HPP

//---------------------------------------------------------------------------//
// end DTK_Exception.hpp
//---------------------------------------------------------------------------//

