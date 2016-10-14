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
 * \brief DTK_ReferenceHexIntegrationRule.cpp
 * \author Stuart R. Slattery
 * \brief Hex reference integration rule.
 */
//---------------------------------------------------------------------------//

#include "DTK_ReferenceHexIntegrationRule.hpp"
#include "DTK_DBC.hpp"

#include <Shards_BasicTopologies.hpp>

namespace DataTransferKit
{
namespace UnitTest
{
//---------------------------------------------------------------------------//
// Constructor.
ReferenceHexIntegrationRule::ReferenceHexIntegrationRule()
    : d_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Given an entity and an integration order, get its integration rule.
void ReferenceHexIntegrationRule::getIntegrationRule(
    const DataTransferKit::Entity &entity, const int order,
    Teuchos::Array<Teuchos::Array<double>> &reference_points,
    Teuchos::Array<double> &weights ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() );
    d_intrepid_rule.getIntegrationRule( d_topo, order, reference_points,
                                        weights );
}

//---------------------------------------------------------------------------//

} // end namespace UnitTest
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ReferenceHexIntegrationRule.hpp
//---------------------------------------------------------------------------//
