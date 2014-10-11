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
 * \brief DTK_MapOperator.hpp
 * \author Stuart R. Slattery
 * \brief Map operator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAPOPERATOR_HPP
#define DTK_MAPOPERATOR_HPP

#include "DTK_FunctionSpace.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_LinearOpBase.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MapOperator
  \brief Map operator interface.

  A map operator maps a field in one entity set to another entity set.
*/
//---------------------------------------------------------------------------//
class MapOperator
{
  public:

    /*!
     * \brief Constructor.
     */
    MapOperator();

    /*!
     * \brief Destructor.
     */
    virtual ~MapOperator();

    //@{
    //! Setup functions.
    /*
     * \brief Setup the map operator from a domain entity set and a range
     * entity set.
     * \param domain_function The field function that contains the data that
     * will be sent to the range.
     * \param range_set The field function that will receive the data from the
     * domain.
     * \param parameters Parameters for the setup.
     */
    virtual void setup( const Teuchos::RCP<FunctionSpace>& domain_space,
			const Teuchos::RCP<FunctionSpace>& range_space,
			const Teuchos::ParameterList& parameters );
    //@}

    //@{
    //! Apply functions.
    /*
     * \brief Apply the map operator to data defined on the entities by
     * computing g = Minv*(v-A*f)
     * \param domain DOFs defined on the domain entities that will be sent to
     * the range.
     * \param range DOFs defined on the range entities that will be
     * received from the domain.
     */
    virtual void apply( 
	const Thyra::MultiVectorBase<double>& domain_dofs,
	const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& range_dofs,
	const double alpha,
	const double beta ) const;
    //@}

  protected:

    //! Mass matrix inverse.
    Teuchos::RCP<Thyra::LinearOpBase<double> > b_mass_matrix_inv;

    //! Coupling matrix.
    Teuchos::RCP<Thyra::LinearOpBase<double> > b_coupling_matrix;

    //! Forcing vector.
    Teuchos::RCP<Thyra::MultiVectorBase<double> > b_forcing_vector;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MAPOPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MapOperator.hpp
//---------------------------------------------------------------------------//
