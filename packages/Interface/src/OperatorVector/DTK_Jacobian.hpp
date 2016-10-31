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
 * \brief DTK_Jacobian.hpp
 * \author Stuart R. Slattery
 * \brief Map operator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_JACOBIAN_HPP
#define DTK_JACOBIAN_HPP

#include <Teuchos_RCP.hpp>

#include "DTK_FunctionSpace.hpp"
#include "DTK_Types.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class Jacobian
  \brief A wrapper implementing Jacobian functionality.
*/
//---------------------------------------------------------------------------//
class Jacobian
{
  public:
    /*!
     * \brief Constructor.
     */
    Jacobian( const Teuchos::RCP<FunctionSpace> &function_space );

    /*!
     * \brief Destructor.
     */
    virtual ~Jacobian();

    /*!
     * \brief Get the Jacobian for a reference point in an entity
     *
     * \param entity The reference point is associated with this entity.
     *
     * \param reference_point The reference point to calculate the Jacobian for.
     *
     * \param[out] jac The computed Jacobian.
     */
    void jacobian( const Entity &entity,
                   const Teuchos::ArrayView<const double> &reference_point,
                   Teuchos::Array<Teuchos::Array<double>> &jac ) const;

    /*!
     * \brief Get the determinant of the Jacobian for a reference point in an
     * entity
     *
     * \param entity The reference point is associated with this entity.
     *
     * \param reference_point The reference point to calculate the Jacobian for.
     *
     * \return The computed determinant.
     */
    double jacobian_determinant(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &reference_point ) const;

  private:
    // The function space for entities.
    Teuchos::RCP<FunctionSpace> d_function_space;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_JACOBIAN_HPP

//---------------------------------------------------------------------------//
// end DTK_Jacobian.hpp
//---------------------------------------------------------------------------//
