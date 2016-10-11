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
 * \brief DTK_MapOperator.cpp
 * \author Stuart R. Slattery
 * \brief Map operator interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_MapOperator.hpp"
#include "DTK_DBC.hpp"
#include "DTK_FieldMultiVector.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MapOperator::MapOperator( const Teuchos::RCP<const TpetraMap> &domain_map,
                          const Teuchos::RCP<const TpetraMap> &range_map )
    : d_domain_map( domain_map )
    , d_range_map( range_map )
    , d_setup_is_complete( false )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Destructor.
MapOperator::~MapOperator() { /* ... */}

//---------------------------------------------------------------------------//
// Get the range map.
Teuchos::RCP<const typename MapOperator::TpetraMap>
MapOperator::getDomainMap() const
{
    DTK_REQUIRE( Teuchos::nonnull( d_domain_map ) );
    return d_domain_map;
}

//---------------------------------------------------------------------------//
// Get the domain map.
Teuchos::RCP<const typename MapOperator::TpetraMap>
MapOperator::getRangeMap() const
{
    DTK_REQUIRE( Teuchos::nonnull( d_range_map ) );
    return d_range_map;
}

//---------------------------------------------------------------------------//
// Setup the map operator.
void MapOperator::setup( const Teuchos::RCP<FunctionSpace> &domain_space,
                         const Teuchos::RCP<FunctionSpace> &range_space )
{
    setupImpl( domain_space, range_space );
    d_setup_is_complete = true;
}

//---------------------------------------------------------------------------//
bool MapOperator::setupIsComplete() const { return d_setup_is_complete; }

//---------------------------------------------------------------------------//
// Apply the map operator.
void MapOperator::apply( const TpetraMultiVector &X, TpetraMultiVector &Y,
                         Teuchos::ETransp mode, const double alpha,
                         const double beta ) const
{
    // Pull data from the applications.
    const FieldMultiVector &X_fmv = dynamic_cast<const FieldMultiVector &>( X );
    const_cast<FieldMultiVector &>( X_fmv ).pullDataFromApplication();
    dynamic_cast<FieldMultiVector &>( Y ).pullDataFromApplication();

    // Apply the operator.
    applyImpl( X, Y, mode, alpha, beta );

    // Push the data into the application.
    dynamic_cast<FieldMultiVector &>( Y ).pushDataToApplication();
}

//---------------------------------------------------------------------------//
// Check if the map has a transpose apply option.n
bool MapOperator::hasTransposeApply() const { return hasTransposeApplyImpl(); }

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MapOperator.cpp
//---------------------------------------------------------------------------//
