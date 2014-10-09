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
 * \brief DTK_MapOperator.cpp
 * \author Stuart R. Slattery
 * \brief Map operator interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_MapOperator.hpp"
#include "DTK_DBC.hpp"

#include <Thyra_MultiVectorStdOps.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MapOperator::MapOperator()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
MapOperator::~MapOperator()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return a string indicating the derived map operator type.
std::string MapOperator::name() const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
    return "Not Implemented";
}

//---------------------------------------------------------------------------//
// Setup the map operator from a source entity set and a target entity set.
void MapOperator::setup( const Teuchos::RCP<FieldGroup>& source_group,
			 const Teuchos::RCP<FieldGroup>& target_group,
			 const Teuchos::ParameterList>& parameters )
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Apply the map operator to data defined on the entities.
void MapOperator::apply( const std::string& source_field_name,
			 const std::string& target_field_name ) const
{
    DTK_REQUIRE( Teuchos::nonnull(b_coupling_matrix) );
    DTK_REQUIRE( Teuchos::nonnull(source_field) );
    DTK_REQUIRE( Teuchos::nonnull(target_field) );

    // Build the source and target vectors.
    Teuchos::RCP<const Thyra::MultiVectorBase<double> >& source_data;
    Teuchos::RCP<Thyra::MultiVectorBase<double> >& target_data;
    buildSourceAndTargetVectors( source_data, target_data );
    DTK_CHECK( Teuchos::nonnull(source_data) );
    DTK_CHECK( Teuchos::nonnull(target_data) );

    // Compute g = Minv*(v-A*f)
    Teuchos::RCP<Thyra::MultiVectorBase<double> > work = 
	b_forcing_vector->clone_mv();
    b_coupling_matrix->apply( 
	Thyra::NOTRANS, *source_data, Teuchos::ptr(work.getRawPtr()), 1.0, 1.0 );
    Thyra::Vt_S( Teuchos::ptr(work.getRawPtr()), -1.0 );
    if ( Teuchos::nonnull(b_forcing_vector) )
    {
	Thyra::Vp_V( Teuchos::ptr(work.getRawPtr()), *b_forcing_vector );
    }
    if ( Teuchos::nonnull(b_mass_matrix_inv) )
    {
	b_mass_matrix_inv->apply( 
	    Thyra::NOTRANS, *work, Teuchos::ptr(target_data.getRawPtr()), 1.0, 1.0 );
    }
    else
    {
	target_data = work;
    }
}

//---------------------------------------------------------------------------//
void MapOperator::buildSourceAndTargetVectors( 
    const Teuchos::RCP<const Field>& source_field,
    const Teuchos::RCP<Field>& target_field,
    Teuchos::RCP<const Thyra::MultiVectorBase<double> >& source_data,
    Teuchos::RCP<Thyra::MultiVectorBase<double> >& target_data ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MapOperator.cpp
//---------------------------------------------------------------------------//
