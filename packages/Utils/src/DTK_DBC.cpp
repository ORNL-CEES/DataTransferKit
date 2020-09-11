/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
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
 * \param file A string containing the file name in which the assertion failed.
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
 * \param file A string containing the file name in which the assertion failed.
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
 * \param file A string containing the file name in which the assertion failed.
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
/*!
 * \brief Throw a DataTransferKit::DataTransferKitException when a user
 * function is missing.
 *
 * \param cond The missing function.
 */
void missingUserFunction( const std::string &cond )
{
    std::ostringstream output_msg;
    output_msg << "Missing user function : " << cond << std::endl;
    throw DataTransferKitException( output_msg.str() );
}

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_DBC.cpp
//---------------------------------------------------------------------------//
