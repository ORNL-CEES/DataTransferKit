//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_MappingStatus.hpp
 * \author Stuart R. Slattery
 * \brief Mapping status.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAPPINGSTATUS_HPP
#define DTK_MAPPINGSTATUS_HPP

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MappingStatus
  \brief Mapping status provides status indications for parametric mapping
  operations.
*/
//---------------------------------------------------------------------------//
class MappingStatus
{
  public:

    /*!
     * \brief Constructor.
     */
    MappingStatus();

    /*!
     * \brief Data constructor.
     */
    MappingStatus( bool success, unsigned number_of_iterations );

    /*!
     * \brief Destructor.
     */
    ~MappingStatus();

    /*!
     * \brief Indicate success.
     */
    void mappingSucceeded();

    /*!
     * \brief Indicate failure.
     */
    void mappingFailed();
    
    /*!
     * \brief Increment the iteration count.
     */
    void incrementIterations();

    /*!
     * \brief Mapping success status.
     * \return True if the mapping succeeded. False if it failed.
     */
    bool success() const;
    
    /*!
     * \brief Iteration count.
     * \return The number of iterations required to construct the mapping.
     */
    int numberOfIterations() const;

    /*!
     * \brief Reset the mapping status.
     */
    void reset();

  private:

    // Success boolean for last mapping.
    bool d_success;

    // Number of iterations for last mapping.
    unsigned d_num_iters;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MAPPINGSTATUS_HPP

//---------------------------------------------------------------------------//
// end DTK_MappingStatus.hpp
//---------------------------------------------------------------------------//
