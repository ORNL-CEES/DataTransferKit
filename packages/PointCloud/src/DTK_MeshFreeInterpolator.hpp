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
 * \file DTK_MeshFreeInterpolator.hpp
 * \author Stuart R. Slattery
 * \brief Interface for mesh free interpolation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHFREEINTERPOLATOR_HPP
#define DTK_MESHFREEINTERPOLATOR_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class MeshFreeInterpolator
 * \brief Base class for mesh free interpolation.
 *
 * The MeshFreeInterpolator is the top-level driver for parallel interpolation
 * problems via compactly supported radial basis functions.
 */
//---------------------------------------------------------------------------//
class MeshFreeInterpolator
{
  public:

    // Constructor.
    MeshFreeInterpolator()
    { /* ... */ }

    //! Destructor.
    virtual ~MeshFreeInterpolator()
    { /* ... */ }

    /*!
     * \brief Set the interpolation problem.
     *
     * \param source_centers Coordinates of the source centers. Must be of the
     * same dimension as the interpolator. Coordinates must be strided such
     * that in 3D: (x1,y1,z1,x2,y2,z2,.....,xN,yN,zN). If there are no source
     * centers on this process then provide an array view of size 0.
     *
     * \param target_centers Coordinates of the target centers. Must be of the
     * same dimension as the interpolator. Coordinates must be strided such
     * that in 3D: (x1,y1,z1,x2,y2,z2,...,xN,yN,zN). If there are no target
     * centers on this process then provide an array view of size 0.
     *
     * \param radius Radius of support for the interpolation. Must be greater
     * than 0.
     */
    virtual void setProblem(
	const Teuchos::ArrayView<const double>& source_centers,
	const Teuchos::ArrayView<const double>& target_centers,
	const double radius ) = 0;

    /*!
     * \brief Given a set of scalar values at the given source centers in the
     * source decomposition, interpolate them onto the target centers in the
     * target decomposition.
     *
     * \param source_data Data values at the source centers. For
     * multi-dimensional data these values must be blocked such that for
     * 3-dimensional data where a, b, and c are the dimensions: 
     * (a1,a2,...,aN,b1,b2,...,bN,c1,c2,...,cN). Length of this array must be
     * dimension * number of source centers. If there is no data on this
     * process then the view must be of size 0.
     *
     * \param target_data Empty data values at the target centers. For
     * multi-dimensional data these values must be blocked such that for
     * 3-dimensional data where a, b, and c are the dimensions:
     * (a1,a2,...,aN,b1,b2,...,bN,c1,c2,...,cN). The interpolated values will
     * be written into this array. Length of this array must be dimension *
     * number of target centers. If there is no data on this process then the
     * view must be of size 0.
     *
     * \param data_dim The dimension of the data. Must be greater than 0.
     */
    virtual void interpolate( 
	const Teuchos::ArrayView<const double>& source_data,
	const Teuchos::ArrayView<double>& target_data,
	const int data_dim ) const = 0;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MESHFREEINTERPOLATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshFreeInterpolator.hpp
//---------------------------------------------------------------------------//

