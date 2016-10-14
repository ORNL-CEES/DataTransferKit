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
 * \brief DTK_ReferenceHexShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief Reference hex shape function
 */
//---------------------------------------------------------------------------//

#include "DTK_ReferenceHexShapeFunction.hpp"
#include "DTK_DBC.hpp"
#include "DTK_ReferenceHexImpl.hpp"
#include "DTK_ReferenceNodeImpl.hpp"

#include <Shards_BasicTopologies.hpp>

namespace DataTransferKit
{
namespace UnitTest
{
//---------------------------------------------------------------------------//
// Constructor.
ReferenceHexShapeFunction::ReferenceHexShapeFunction()
    : d_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Given an entity, get the ids of the degrees of freedom in the vector space
// supporting its shape function.
void ReferenceHexShapeFunction::entitySupportIds(
    const DataTransferKit::Entity &entity,
    Teuchos::Array<DataTransferKit::SupportId> &support_ids ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() ||
                 0 == entity.topologicalDimension() );

    // Node case.
    if ( 0 == entity.topologicalDimension() )
    {
        support_ids.resize( 1 );
        support_ids[0] = Teuchos::rcp_dynamic_cast<ReferenceNodeExtraData>(
                             entity.extraData() )
                             ->id;
    }

    // Hex case.
    else
    {
        support_ids.resize( 8 );
        auto &node_ids = Teuchos::rcp_dynamic_cast<ReferenceHexExtraData>(
                             entity.extraData() )
                             ->node_ids;
        DTK_CHECK( 8 == node_ids.size() );
        std::copy( node_ids.begin(), node_ids.end(), support_ids.begin() );
    }
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void ReferenceHexShapeFunction::evaluateValue(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<double> &values ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() );
    d_intrepid_shape.evaluateValue( d_topo, reference_point, values );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void ReferenceHexShapeFunction::evaluateGradient(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<Teuchos::Array<double>> &gradients ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() );
    d_intrepid_shape.evaluateGradient( d_topo, reference_point, gradients );
}

//---------------------------------------------------------------------------//

} // end namespace UnitTest
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ReferenceHexShapeFunction.cpp
//---------------------------------------------------------------------------//
