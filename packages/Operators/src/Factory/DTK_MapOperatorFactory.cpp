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
 * \brief DTK_MapOperatorFactory.cpp
 * \author Stuart R. Slattery
 * \brief Map operator factory.
 */
//---------------------------------------------------------------------------//

#include "DTK_MapOperatorFactory.hpp"
#include "DTK_DBC.hpp"
#include "DTK_L2ProjectionOperator.hpp"
#include "DTK_ConsistentInterpolationOperator.hpp"
#include "DTK_PointCloudOperatorFactory.hpp"

namespace DataTransferKit
{
//--------------------------------------------------------------------------//
// Constructor
MapOperatorFactory::MapOperatorFactory()
{
    d_name_map["L2 Projection"] = L2_PROJECTION;
    d_name_map["Consistent Interpolation"] = CONSISTENT_INTERPOLATION;
    d_name_map["Point Cloud"] = POINT_CLOUD;
}
    
//---------------------------------------------------------------------------//
// Creation method.
Teuchos::RCP<MapOperator>
MapOperatorFactory::create(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map,
    const Teuchos::ParameterList& parameters )
{
    // Get the name of the map to build.
    std::string map_name = parameters.get<std::string>("Map Type");
    DTK_REQUIRE( d_name_map.count(map_name) );
    int map_id = d_name_map.find( map_name )->second;

    // Get the parameters for the map.
    Teuchos::ParameterList map_list = parameters.sublist(map_name);

    // Initialize subclass factories.
    PointCloudOperatorFactory pcloud_factory;

    // Build the map.
    Teuchos::RCP<MapOperator> map;
    switch( map_id )
    {
	// L2 Projection.
	case L2_PROJECTION:
	    map = Teuchos::rcp(	
		new L2ProjectionOperator(domain_map,range_map,parameters) );
	    break;
	
	// Consistent Interpolation.
	case CONSISTENT_INTERPOLATION:
	    map = Teuchos::rcp(	new ConsistentInterpolationOperator(
				    domain_map,range_map,parameters) );
	    break;

	// Spline Interpolation
	case POINT_CLOUD:
	    map = pcloud_factory.create( domain_map, range_map, map_list );
	    break;

        // Throw on default.
	default:
	    bool map_type_is_valid = false;
	    DTK_INSIST( map_type_is_valid );
	    break;
    }

    DTK_ENSURE( Teuchos::nonnull(map) );
    return map;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MapOperatorFactory.cpp
//---------------------------------------------------------------------------//
