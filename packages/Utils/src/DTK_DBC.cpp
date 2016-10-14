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
 * \file   DTK_DBC.cpp
 * \author Stuart Slattery
 * \brief  Assertions for error handling and Design-by-Contract.
 */
//---------------------------------------------------------------------------//

#include <sstream>

#include "DTK_DBC.hpp"

#include <Teuchos_stacktrace.hpp>

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
 * \return DataTransferKitException output.
 */
std::string DataTransferKitException::generate_output( const std::string &cond,
                                                       const std::string &file,
                                                       const int line ) const
{
    std::ostringstream output;
    output << "DataTransferKit DataTransferKitException: " << cond
           << ", failed in " << file << ", line " << line << "." << std::endl;
    return output.str();
}

//---------------------------------------------------------------------------//
// Throw functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Throw a DataTransferKit::DataTransferKitException.
 *
 * \param cond A string containing the assertion condition that failed.
 *
 * \param field A string containing the file name in which the assertion
 * failed.
 *
 * \param line The line number at which the assertion failed.
 */
void throwDataTransferKitException( const std::string &cond,
                                    const std::string &file, const int line )
{
#ifdef HAVE_TEUCHOS_STACKTRACE
    // If Teuchos stacktrace is turned on, store the stack before we throw so
    // we can get it later.
    Teuchos::store_stacktrace();
#endif
    throw DataTransferKitException( cond, file, line );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Throw a DataTransferKit::DataTransferKitException when an error code
 * fails.
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
void errorCodeFailure( const std::string &cond, const std::string &file,
                       const int line, const int error_code )
{
#ifdef HAVE_TEUCHOS_STACKTRACE
    // If Teuchos stacktrace is turned on, store the stack before we throw so
    // we can get it later.
    Teuchos::store_stacktrace();
#endif
    std::ostringstream output_msg;
    output_msg << "Error code : " << cond << ", failed in " << file << ":"
               << line << std::endl
               << "with error code:" << std::endl
               << "\"" << error_code << "\"" << std::endl;
    throw DataTransferKitException( output_msg.str() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_DBC.cpp
//---------------------------------------------------------------------------//
