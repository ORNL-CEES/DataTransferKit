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
 * \file DTK_BoundingVolumeList.hpp
 * \brief BoundingVolume list.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOUNDINGVOLUMELIST_HPP
#define DTK_BOUNDINGVOLUMELIST_HPP

#include "DTK_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class BoundingVolumeList.
 *
 * \brief Trivially-copyable bounding volume list.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 */
template <class... ViewProperties>
class BoundingVolumeList
{
  public:
    //! View traits.
    using ViewTraits = typename Kokkos::ViewTraits<int, ViewProperties...>;

    //! The coordinates of the bounding volumes that are locally-owned by this
    //! MPI rank. This view is rank-3 and should be sized as (number of
    //! bounding volumes, spatial dimension, 2 [for low and high coordinate]).
    Kokkos::View<Coordinate * * [2], ViewProperties...> bounding_volumes;

    //! View indicating if the given bounding volume is owned by the local
    //! process or is a ghost. This information is not necessary for all
    //! algorithms and therefore this view can optionally be of size 0 or NULL
    //! to indicate that no ghost information is available thereby indicating
    //! that all bounding volumes provided are locally-owned by this MPI
    //! rank. This view is rank-1 and of length of the number of bounding
    //! volumes in the list.
    Kokkos::View<bool *, ViewProperties...> is_ghost_volume;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_BOUNDINGVOLUMELIST_HPP

//---------------------------------------------------------------------------//
// end DTK_BoundingVolumeList.hpp
//---------------------------------------------------------------------------//
