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
 * \file DTK_MeshFreeInterpolatorFactory_impl.hpp
 * \author Stuart R. Slattery
 * \brief Linear solver factory implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHFREEINTERPOLATORFACTORY_IMPL_HPP
#define DTK_MESHFREEINTERPOLATORFACTORY_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_SplineInterpolator.hpp"
#include "DTK_MovingLeastSquare.hpp"
#include "DTK_WendlandBasis.hpp"
#include "DTK_WuBasis.hpp"
#include "DTK_BuhmannBasis.hpp"

namespace DataTransferKit
{


//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
MeshFreeInterpolatorFactory::MeshFreeInterpolatorFactory()
{
    // Create the interpolator name-to-enum map.
    d_interpolator_map["Spline"] = SPLINE;
    d_interpolator_map["Moving Least Square"] = MOVING_LEAST_SQUARE;

    // Create the basis name-to-enum map.
    d_basis_map["Wendland"] = WENDLAND;
    d_basis_map["Wu"] = WU;
    d_basis_map["Buhmann"] = BUHMANN;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Creation method.
 */
template<class GO>
Teuchos::RCP<MeshFreeInterpolator>
MeshFreeInterpolatorFactory::create( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const std::string& interplator_name, 
    const std::string& basis_name,
    const int basis_order,
    const int space_dim )
{
    DTK_REQUIRE( Teuchos::nonnull(comm) );

    Teuchos::RCP<MeshFreeInterpolator> interpolator;

    MapType::const_iterator interpolator_id = 
	d_interpolator_map.find( interpolator_name );
    DTK_INSIST( id != d_interpolator_map.end() );

    MapType::const_iterator basis_id = d_basis_map.find( basis_name );
    DTK_INSIST( id != d_basis_map.end() );

    switch( interpolator_id->second )
    {
	case SPLINE:
	{
	    switch( basis_id->second )
	    {
		case WENDLAND:
		{
		    switch basis_order:
		    {
			case 0:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<0>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<0>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<0>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			case 2:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<2>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<2>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<2>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;
			
			case 4:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<4>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<4>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<4>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			case 6:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<6>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<6>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WendlandBasis<6>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			default:
			    throw Assertion("Basis order not supported!");
			    break;
		    }
		}
		break;

		case WU:
		{
			case 2:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WuBasis<2>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WuBasis<2>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WuBasis<2>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;
			
			case 4:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WuBasis<4>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WuBasis<4>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					WuBasis<4>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			default:
			    throw Assertion("Basis order not supported!");
			    break;
		}
		break;
		
		case BUHMANN:
		{
			case 3:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					BuhmannBasis<3>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					BuhmannBasis<3>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new SplineInterpolator<
					BuhmannBasis<3>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			default:
			    throw Assertion("Basis order not supported!");
			    break;
		}
		break;

		default:
		    throw Assertion("Basis type not supported!");
		    break;
	    }
	}
	break;

	case MOVING_LEAST_SQUARE:
	{
	    switch( basis_id->second )
	    {
		case WENDLAND:
		{
		    switch basis_order:
		    {
			case 0:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<0>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<0>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<0>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			case 2:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<2>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<2>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<2>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;
			
			case 4:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<4>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<4>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<4>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			case 6:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<6>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<6>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WendlandBasis<6>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			default:
			    throw Assertion("Basis order not supported!");
			    break;
		    }
		}
		break;

		case WU:
		{
			case 2:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WuBasis<2>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WuBasis<2>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WuBasis<2>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;
			
			case 4:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WuBasis<4>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WuBasis<4>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					WuBasis<4>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			default:
			    throw Assertion("Basis order not supported!");
			    break;
		}
		break;
		
		case BUHMANN:
		{
			case 3:
			{
			    switch space_dim:
			    {
				case 1:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					BuhmannBasis<3>,GO,1>(comm) );
				}
				break;
				case 2:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					BuhmannBasis<3>,GO,2>(comm) );
				}
				break;
				case 3:
				{
				    interpolator = Teuchos::rcp(
					new MovingLeastSquare<
					BuhmannBasis<3>,GO,3>(comm) );
				}
				break;

				default:
				    throw Assertion("Spatial Dimension not supported!");
				    break;
			    }
			}
			break;

			default:
			    throw Assertion("Basis order not supported!");
			    break;
		}
		break;

		default:
		    throw Assertion("Basis type not supported!");
		    break;
	    }
	}
	break;

	default:
	    throw Assertion("Interpolator type not supported!");
	    break;
    }

    DTK_ENSURE( Teuchos::nonnull(interpolator) );
    return interpolator;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESHFREEINTERPOLATORFACTORY_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshFreeInterpolatorFactory_impl.hpp
// ---------------------------------------------------------------------------//

