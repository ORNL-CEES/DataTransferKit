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
 * \brief Mesh free interpolator factory implementation.
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
 * \brief Creation method.
 *
 * \param comm The parallel communicator over which to build the interpolator.
 *
 * \param interpolator_name The name of the interpolation type.
 *
 * \param basis_name The name of the basis type supporting the interpolation.
 *
 * \param basis_order The order of the basis supporting the interpolation.
 *
 * \param space_dim The spatial dimension of the interpolation problem.
 */
template<class GO>
Teuchos::RCP<MeshFreeInterpolator>
MeshFreeInterpolatorFactory::create( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const std::string& interpolator_name, 
    const std::string& basis_name,
    const int basis_order,
    const int space_dim )
{
    DTK_REQUIRE( Teuchos::nonnull(comm) );

    Teuchos::RCP<MeshFreeInterpolator> interpolator;

    int interpolator_id = -1;
    if ( "Spline" == interpolator_name ) interpolator_id = 0;
    else if ( "Moving Least Square" == interpolator_name ) interpolator_id = 1;

    int basis_id = -1;
    if ( "Wendland" == basis_name ) basis_id = 0;
    else if ( "Wu" == basis_name ) basis_id = 1;
    else if ( "Buhmann" == basis_name ) basis_id = 2;

    switch( interpolator_id )
    {
	// Spline
	case 0:
	{
	    // Wendland
	    switch( basis_id )
	    {
		case 0:
		{
		    switch( basis_order )
		    {
			case 0:
			{
			    switch( space_dim )
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
			    switch( space_dim )
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
			    switch( space_dim )
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
			    switch( space_dim )
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

		// Wu
		case 1:
		{
		    switch( basis_order )
		    {
			case 2:
			{
			    switch( space_dim )
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
			    switch( space_dim )
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
		}
		break;
		
		// Buhmann
		case 2:
		{
		    switch( basis_order )
		    {
			case 3:
			{
			    switch( space_dim )
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
		}
		break;

		default:
		    throw Assertion("Basis type not supported!");
		    break;
	    }
	}
	break;

	// Moving Least Square.
	case 1:
	{
	    switch( basis_id )
	    {
		// Wendland
		case 0:
		{
		    switch( basis_order )
		    {
			case 0:
			{
			    switch( space_dim )
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
			    switch( space_dim )
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
			    switch( space_dim )
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
			    switch( space_dim )
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

		// Wu
		case 1:
		{
		    switch( basis_order )
		    {
			case 2:
			{
			    switch( space_dim )
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
			    switch( space_dim )
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
		}
		break;
		
		// Buhmann
		case 2:
		{
		    switch( basis_order )
		    {
			case 3:
			{
			    switch( space_dim )
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

