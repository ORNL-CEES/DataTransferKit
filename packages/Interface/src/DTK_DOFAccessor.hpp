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
 * \brief DTK_DOFAccessor.hpp
 * \author Stuart R. Slattery
 * \brief Interface for accessing DOF data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DOFACCESSOR_HPP
#define DTK_DOFACCESSOR_HPP

#include <string>

#include "DTK_Entity.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class DOFAccessor
  \brief Interface for accessing DOF data for a field.
*/
//---------------------------------------------------------------------------//
class DOFAccessor
{
  public:

    /*!
     * \brief Constructor.
     */
    DOFAccessor();

    /*!
     * \brief Destructor.
     */
    virtual ~DOFAccessor();
    
    /*!
     * \brief Get the number of dimensions in the field.
     * \return The dimension of the field.
     */
    virtual int dimension() const;

    /*!
     * \brief Get the size of each dimension in the field.
     * \param dims The size of each dimension in the field.
     */
    virtual void dimensionSizes( Teuchos::Array<int>& dim_sizes ) const;

    /*!
     * \brief Given a degree of freedom id, get a const view of the data for
     * that degree of freedom. 
     * \param dof_id Get the data for the degree of freedom with this id.
     * \param dof_values A const view of the degree of freedom data.
     */
    virtual void access( const EntityId& dof_id,
			 Kokkos::View<const double>& dof_values ) const;

    /*!
     * \brief Given a degree of freedom id, get a nonconst view of the data for
     * that degree of freedom. 
     * \param dof_id Get the data for the degree of freedom with this id.
     * \param dof_values A nonconst view of the degree of freedom data.
     */
    virtual void accessNonConst( const EntityId& dof_id,
				 Kokkos::View<double>& dof_values ) const;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_DOFACCESSOR_HPP

//---------------------------------------------------------------------------//
// end DTK_DOFAccessor.hpp
//---------------------------------------------------------------------------//
