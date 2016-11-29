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
 * \brief DTK_ReferenceHexLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include "DTK_ReferenceHexLocalMap.hpp"
#include "DTK_DBC.hpp"
#include "DTK_IntrepidCellLocalMap.hpp"
#include "DTK_ReferenceHexImpl.hpp"
#include "DTK_ReferenceNodeImpl.hpp"

#include <Shards_BasicTopologies.hpp>

namespace DataTransferKit
{
namespace UnitTest
{
//---------------------------------------------------------------------------//
// Constructor.
ReferenceHexLocalMap::ReferenceHexLocalMap()
    : d_inclusion_tol( 1.0e-6 )
    , d_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Set parameters for mapping.
void ReferenceHexLocalMap::setParameters(
    const Teuchos::ParameterList &parameters )
{
    if ( parameters.isParameter( "Point Inclusion Tolerance" ) )
    {
        d_inclusion_tol = parameters.get<double>( "Point Inclusion Tolerance" );
    }
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
double
ReferenceHexLocalMap::measure( const DataTransferKit::Entity &entity ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() ||
                 0 == entity.topologicalDimension() );

    // Node case.
    if ( 0 == entity.topologicalDimension() )
    {
        return 0.0;
    }

    // Hex case.
    else
    {
        auto &cell_coords = Teuchos::rcp_dynamic_cast<ReferenceHexExtraData>(
                                entity.extraData() )
                                ->node_coords;
        return DataTransferKit::IntrepidCellLocalMap::measure( d_topo,
                                                               cell_coords );
    }
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void ReferenceHexLocalMap::centroid(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<double> &centroid ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() ||
                 0 == entity.topologicalDimension() );

    // Node case.
    if ( 0 == entity.topologicalDimension() )
    {
        auto &node_coords = Teuchos::rcp_dynamic_cast<ReferenceNodeExtraData>(
                                entity.extraData() )
                                ->node_coords;
        std::copy( node_coords.begin(), node_coords.end(), centroid.begin() );
    }

    // Hex case.
    else
    {
        auto &cell_coords = Teuchos::rcp_dynamic_cast<ReferenceHexExtraData>(
                                entity.extraData() )
                                ->node_coords;
        DataTransferKit::IntrepidCellLocalMap::centroid( d_topo, cell_coords,
                                                         centroid );
    }
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
bool ReferenceHexLocalMap::mapToReferenceFrame(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &physical_point,
    const Teuchos::ArrayView<double> &reference_point ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() );

    auto &cell_coords =
        Teuchos::rcp_dynamic_cast<ReferenceHexExtraData>( entity.extraData() )
            ->node_coords;
    return DataTransferKit::IntrepidCellLocalMap::mapToReferenceFrame(
        d_topo, cell_coords, physical_point, reference_point );
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool ReferenceHexLocalMap::checkPointInclusion(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() );

    return DataTransferKit::IntrepidCellLocalMap::checkPointInclusion(
        d_topo, reference_point, d_inclusion_tol );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void ReferenceHexLocalMap::mapToPhysicalFrame(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &physical_point ) const
{
    DTK_REQUIRE( 3 == entity.topologicalDimension() );

    auto &cell_coords =
        Teuchos::rcp_dynamic_cast<ReferenceHexExtraData>( entity.extraData() )
            ->node_coords;
    DataTransferKit::IntrepidCellLocalMap::mapToPhysicalFrame(
        d_topo, cell_coords, reference_point, physical_point );
}

//---------------------------------------------------------------------------//

} // end namespace UnitTest
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ReferenceHexLocalMap.cpp
//---------------------------------------------------------------------------//
