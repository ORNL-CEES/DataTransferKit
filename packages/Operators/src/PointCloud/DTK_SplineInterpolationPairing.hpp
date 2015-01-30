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
 * \file   DTK_SplineInterpolationPairing.hpp
 * \author Stuart R. Slattery
 * \brief  Local source/parent center pairings.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTERPOLATIONPAIRING_HPP
#define DTK_INTERPOLATIONPAIRING_HPP

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class SplineInterpolationPairing
 * \brief Local child/parent center pairings.
 *
 * Build groups of local child centers that are within the given radius of
 * the parent centers. Each parent center will have a list of child centers.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class SplineInterpolationPairing
{
  public:

    // Constructor.
    SplineInterpolationPairing( 
	const Teuchos::ArrayView<const double>& child_centers,
	const Teuchos::ArrayView<const double>& parent_centers,
	const double radius );

    // Given a parent center local id get the ids of the child centers within
    // the given radius.
    Teuchos::ArrayView<const unsigned> 
    childCenterIds( const unsigned parent_id ) const;

    // Get the number of child centers per parent center.
    Teuchos::ArrayRCP<std::size_t> childrenPerParent() const
    { return d_pair_sizes; }

  private:

    // Pairings.
    Teuchos::Array<Teuchos::Array<unsigned> > d_pairings;

    // Number of child centers per parent center.
    Teuchos::ArrayRCP<std::size_t> d_pair_sizes;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_SplineInterpolationPairing_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_INTERPOLATIONPAIRING_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineInterpolationPairing.hpp
//---------------------------------------------------------------------------//

