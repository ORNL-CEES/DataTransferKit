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
 * \brief DTK_MapOperatorFactory_impl.hpp
 * \author Stuart R. Slattery
 * \brief Map operator factory.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAPOPERATORFACTORY_IMPL_HPP
#define DTK_MAPOPERATORFACTORY_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_ConsistentInterpolationOperator.hpp"
#include "DTK_PointCloudOperatorFactory.hpp"

namespace DataTransferKit
{
//--------------------------------------------------------------------------//
// Constructor
template<class Scalar>
MapOperatorFactory<Scalar>::MapOperatorFactory()
{
    d_name_map["Consistent Interpolation"] = CONSISTENT_INTERPOLATION;
    d_name_map["Point Cloud"] = POINT_CLOUD;
}
    
//---------------------------------------------------------------------------//
// Creation method.
template<class Scalar>
Teuchos::RCP<MapOperator<Scalar> >
MapOperatorFactory::create(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map,
    const Teuchos::ParameterList& parameters )
{
    DTK_REQUIRE( Teuchos::nonnull(parameters) );

    // Get the name of the map to build.
    std::string map_name = parameters.get<std::string>("Map Type");
    DTK_REQUIRE( d_name_map.count(map_name) );
    int map_id = d_name_map.find( map_name )->second;

    // Get the parameters for the map.
    Teuchos::ParameterList map_list = parameters.sublist(map_name);

    // Build the map.
    Teuchos::RCP<MapOperator<Scalar> > map;
    switch( map_id )
    {
	// Consistent Interpolation.
	case CONSISTENT_INTERPOLATION:
	    map = Teuchos::rcp(	new ConsistentInterpolationOperator<Scalar>(
				    domain_map,range_map) );
	    break;

	// Spline Interpolation
	case POINT_CLOUD:
	    PointCloudOperatorFactory<Scalar> pcloud_factory;
	    map = pcloud_factory.create( domain_map, range_map, map_list );
	    break;
	    
	default:
	    DTK_INSIST( false, "Map operator type not supported!" );
	    break;
    }

    return map;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MAPOPERATORFACTORY_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MapOperator_impl.hpp
//---------------------------------------------------------------------------//
