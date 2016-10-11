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
 * \file   DTK_Assertion.hpp
 * \author Stuart Slattery
 * \brief  Assertions and Design-by-Contract for error handling.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_Assertion_HPP
#define DTK_Assertion_HPP

#include <stdexcept>
#include <string>

#include "DataTransferKitUtils_config.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Base class for DTK assertions. This structure is heavily based on
 * that in Nemesis developed by Tom Evans. We derive from std::logic_error
 * here as the DBC checks that utilize this class are meant to find errors
 * that can be prevented before runtime.
 */
//---------------------------------------------------------------------------//
class Assertion : public std::logic_error
{
  public:
    /*!
     * \brief Default constructor.
     *
     * \param msg Error message.
     */
    Assertion( const std::string &msg )
        : std::logic_error( msg )
    { /* ... */
    }

    /*!
     * \brief DBC constructor.
     *
     * \param cond A string containing the assertion condition that failed.
     *
     * \param field A string containing the file name in which the assertion
     * failed.
     *
     * \param line The line number at which the assertion failed.
     */
    Assertion( const std::string &cond, const std::string &file,
               const int line )
        : std::logic_error( generate_output( cond, file, line ) )
    { /* ... */
    }

    //! Destructor.
    virtual ~Assertion() throw() { /* ... */}

  private:
    // Build an assertion output from advanced constructor arguments.
    std::string generate_output( const std::string &cond,
                                 const std::string &file,
                                 const int line ) const;
};

//---------------------------------------------------------------------------//
// Throw functions.
//---------------------------------------------------------------------------//
// Throw a DataTransferKit::Assertion.
void throwAssertion( const std::string &cond, const std::string &file,
                     const int line );

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Design-by-Contract macros.
//---------------------------------------------------------------------------//
/*!
  \page DataTransferKit Design-by-Contract.

  Design-by-Contract (DBC) functionality is provided to verify function
  preconditions, postconditions, and invariants. These checks are separated
  from the debug build and can be activated for both release and debug
  builds. They can be activated by setting the following in a CMake
  configure:

  -D DataTransferKit_ENABLE_DBC:BOOL=ON

  By default, DBC is deactivated. Although they will require additional
  computational overhead, these checks provide a mechanism for veryifing
  library input arguments. Note that the bounds-checking functionality used
  within the DataTransferKit is only provided by a debug build.

  In addition, DTK_REMEMBER is provided to store values used only for DBC
  checks and no other place in executed code.

  Separate from the DBC build, DTK_INSIST can be used at any time verify a
  conditional. This should be used instead of the standard cassert.

  DTK_CHECK_ERROR_CODE provides DBC support for libraries that return error
  codes with 0 as the value for no errors.
 */

#if HAVE_DTK_DBC

#define testPrecondition( c )                                                  \
    if ( !( c ) )                                                              \
    DataTransferKit::throwAssertion( #c, __FILE__, __LINE__ )
#define testPostcondition( c )                                                 \
    if ( !( c ) )                                                              \
    DataTransferKit::throwAssertion( #c, __FILE__, __LINE__ )
#define testInvariant( c )                                                     \
    if ( !( c ) )                                                              \
    DataTransferKit::throwAssertion( #c, __FILE__, __LINE__ )
#define remember( c ) c

#else

#define testPrecondition( c )
#define testPostcondition( c )
#define testInvariant( c )
#define remember( c )
#endif

#define insist( c )                                                            \
    if ( !( c ) )                                                              \
    DataTransferKit::throwAssertion( #c, __FILE__, __LINE__ )

//---------------------------------------------------------------------------//

#endif // end DTK_Assertion_HPP

//---------------------------------------------------------------------------//
// end DTK_Assertion.hpp
//---------------------------------------------------------------------------//
