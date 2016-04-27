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
 * \brief DTK_PointCloudOperatorFactory.cpp
 * \author Stuart R. Slattery
 * \brief Point cloud map operator factory.
 */
//---------------------------------------------------------------------------//

#include "DTK_PointCloudOperatorFactory.hpp"
#include "DTK_DBC.hpp"
#include "DTK_NodeToNodeOperator.hpp"
#include "DTK_SplineInterpolationOperator.hpp"
#include "DTK_MovingLeastSquareReconstructionOperator.hpp"
#include "DTK_WendlandBasis.hpp"
#include "DTK_WuBasis.hpp"
#include "DTK_BuhmannBasis.hpp"

namespace DataTransferKit
{
//--------------------------------------------------------------------------//
// Constructor
PointCloudOperatorFactory::PointCloudOperatorFactory()
{
    d_name_map["Spline Interpolation"] = SPLINE_INTERPOLATION;
    d_name_map["Moving Least Square Reconstruction"] = MOVING_LEAST_SQUARE;
    d_name_map["Node To Node"] = NODE_TO_NODE;

    d_basis_map["Wendland"] = WENDLAND;
    d_basis_map["Wu"] = WU;
    d_basis_map["Buhmann"] = BUHMANN;
}
    
//---------------------------------------------------------------------------//
// Creation method.
Teuchos::RCP<MapOperator>
PointCloudOperatorFactory::create(
    const Teuchos::RCP<const TpetraMap>& domain_map,
    const Teuchos::RCP<const TpetraMap>& range_map,
    const Teuchos::ParameterList& parameters )
{
    // Get the name of the map to build.
    std::string map_name = parameters.get<std::string>("Map Type");
    DTK_REQUIRE( d_name_map.count(map_name) );
    int map_id = d_name_map.find( map_name )->second;

    // Get the spatial dimension.
    int space_dim = parameters.get<int>("Spatial Dimension");

    // Get the basis type.
    std::string basis_name = parameters.get<std::string>("Basis Type");
    DTK_REQUIRE( d_basis_map.count(basis_name) );
    int basis_id = d_basis_map.find( basis_name )->second;

    // Get the basis order.
    int basis_order = parameters.get<int>("Basis Order");
    
    // Build the map.
    Teuchos::RCP<MapOperator> map;
    switch( map_id )
    {
	// Spline Interpolation.
	case NODE_TO_NODE:
	    switch( space_dim )
	    {
		case 1:
                    map = Teuchos::rcp(
                        new NodeToNodeOperator<1>(
                            domain_map,range_map,parameters) );
                    break;
		case 2:
                    map = Teuchos::rcp(
                        new NodeToNodeOperator<2>(
                            domain_map,range_map,parameters) );
                    break;
		case 3:
                    map = Teuchos::rcp(
                        new NodeToNodeOperator<3>(
                            domain_map,range_map,parameters) );
                    break;
            }
            break;
            
        // Spline Interpolation.
	case SPLINE_INTERPOLATION:
	    switch( space_dim )
	    {
		case 1:
		    switch( basis_id )
		    {
			case WENDLAND:
			    switch( basis_order )
			    {
				case 0:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<0>,1>(domain_map,
                                                            range_map,
                                                            parameters) );
				    break;
				case 2:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<2>,1>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<4>,1>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 6:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<6>,1>(domain_map,
								   range_map,
								   parameters) );
				    break;
			    }
			    break;
			case WU:
			    switch( basis_order )
			    {
				case 2:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WuBasis<2>,1>(domain_map,
							     range_map,
							     parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WuBasis<4>,1>(domain_map,
							     range_map,
							     parameters) );
				    break;
			    }
			    break;
			case BUHMANN:
			    switch( basis_order )
			    {
				case 3:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					BuhmannBasis<3>,1>(domain_map,
								  range_map,
								  parameters) );
				    break;
			    }
			    break;
		    }
		    break;

		case 2:
		    switch( basis_id )
		    {
			case WENDLAND:
			    switch( basis_order )
			    {
				case 0:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<0>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 2:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<2>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<4>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 6:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<6>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
			    }
			    break;
			case WU:
			    switch( basis_order )
			    {
				case 2:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WuBasis<2>,2>(domain_map,
							     range_map,
							     parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WuBasis<4>,2>(domain_map,
							     range_map,
							     parameters) );
				    break;
			    }
			    break;
			case BUHMANN:
			    switch( basis_order )
			    {
				case 3:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					BuhmannBasis<3>,2>(domain_map,
								  range_map,
								  parameters) );
				    break;
			    }
			    break;
		    }
		    break;


		case 3:
		    switch( basis_id )
		    {
			case WENDLAND:
			    switch( basis_order )
			    {
				case 0:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<0>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 2:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<2>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<4>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 6:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WendlandBasis<6>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
			    }
			    break;
			case WU:
			    switch( basis_order )
			    {
				case 2:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WuBasis<2>,3>(domain_map,
							     range_map,
							     parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					WuBasis<4>,3>(domain_map,
							     range_map,
							     parameters) );
				    break;
			    }
			    break;
			case BUHMANN:
			    switch( basis_order )
			    {
				case 3:
				    map = Teuchos::rcp(
					new SplineInterpolationOperator<
					BuhmannBasis<3>,3>(domain_map,
								  range_map,
								  parameters) );
				    break;
			    }
			    break;
		    }
		    break;

	    }
	    break;

	// Moving least square reconstruction.
	case MOVING_LEAST_SQUARE:
	    switch( space_dim )
	    {
		case 1:
		    switch( basis_id )
		    {
			case WENDLAND:
			    switch( basis_order )
			    {
				case 0:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<0>,1>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 2:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<2>,1>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<4>,1>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 6:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<6>,1>(domain_map,
								   range_map,
								   parameters) );
				    break;
			    }
			    break;
			case WU:
			    switch( basis_order )
			    {
				case 2:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WuBasis<2>,1>(domain_map,
							     range_map,
							     parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WuBasis<4>,1>(domain_map,
							     range_map,
							     parameters) );
				    break;
			    }
			    break;
			case BUHMANN:
			    switch( basis_order )
			    {
				case 3:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					BuhmannBasis<3>,1>(domain_map,
								  range_map,
								  parameters) );
				    break;
			    }
			    break;
		    }
		    break;

		case 2:
		    switch( basis_id )
		    {
			case WENDLAND:
			    switch( basis_order )
			    {
				case 0:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<0>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 2:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<2>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<4>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 6:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<6>,2>(domain_map,
								   range_map,
								   parameters) );
				    break;
			    }
			    break;
			case WU:
			    switch( basis_order )
			    {
				case 2:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WuBasis<2>,2>(domain_map,
							     range_map,
							     parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WuBasis<4>,2>(domain_map,
							     range_map,
							     parameters) );
				    break;
			    }
			    break;
			case BUHMANN:
			    switch( basis_order )
			    {
				case 3:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					BuhmannBasis<3>,2>(domain_map,
								  range_map,
								  parameters) );
				    break;
			    }
			    break;
		    }
		    break;


		case 3:
		    switch( basis_id )
		    {
			case WENDLAND:
			    switch( basis_order )
			    {
				case 0:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<0>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 2:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<2>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<4>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
				case 6:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WendlandBasis<6>,3>(domain_map,
								   range_map,
								   parameters) );
				    break;
			    }
			    break;
			case WU:
			    switch( basis_order )
			    {
				case 2:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WuBasis<2>,3>(domain_map,
							     range_map,
							     parameters) );
				    break;
				case 4:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					WuBasis<4>,3>(domain_map,
							     range_map,
							     parameters) );
				    break;
			    }
			    break;
			case BUHMANN:
			    switch( basis_order )
			    {
				case 3:
				    map = Teuchos::rcp(
					new MovingLeastSquareReconstructionOperator<
					BuhmannBasis<3>,3>(domain_map,
								  range_map,
								  parameters) );
				    break;
			    }
			    break;
		    }
		    break;

	    }
	    break;
	    
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
// end DTK_MapOperator.cpp
//---------------------------------------------------------------------------//
