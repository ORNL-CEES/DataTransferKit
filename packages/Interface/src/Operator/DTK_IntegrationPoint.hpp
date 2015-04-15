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
 * \brief DTK_IntegrationPointSet.hpp
 * \author Stuart R. Slattery
 * \brief Integration point set.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTEGRATIONPOINT_HPP
#define DTK_INTEGRATIONPOINT_HPP

#include "DTK_Types.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// IntegrationPoint 
//---------------------------------------------------------------------------//
struct IntegrationPoint
{
    // Id of the entity that owns the point.
    EntityId d_owner_id;

    // Measure of the entity that owns the point.
    double d_owner_measure;
    
    // Weight of the point in the integration rule.
    double d_integration_weight;

    // Physical coordinates of the point in the owning entity.
    Teuchos::Array<double> d_physical_coordinates;

    // Support ids of the owning entity.
    Teuchos::Array<SupportId> d_owner_support_ids;

    // Shape function evaluation of the point in the owning entity.
    Teuchos::Array<double> d_owner_shape_evals;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTEGRATIONPOINT_HPP

//---------------------------------------------------------------------------//
// end DTK_IntegrationPointSet.hpp
//---------------------------------------------------------------------------//
