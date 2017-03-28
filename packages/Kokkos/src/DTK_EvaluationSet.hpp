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
 * \file DTK_EvaluationSet.hpp
 * \brief Evaluation point set.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EVALUATIONSET_HPP
#define DTK_EVALUATIONSET_HPP

#include "DTK_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class EvaluationSet.
 *
 * \brief Trivially-copyable list of evaluation points.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 */
template <class... ViewProperties>
class EvaluationSet
{
  public:
    //@{
    //! Type aliases.

    using ViewTraits = Kokkos::ViewTraits<int, ViewProperties...>;
    using ExecutionSpace = typename ViewTraits::execution_space;
    using MemorySpace = typename ViewTraits::memory_space;
    using Device = typename ViewTraits::device_type;
    //@}

    //! The coordinates of the evaluation points. This view is rank-2 and
    //! should be sized as (number of eval points, spatial dimension)
    Kokkos::View<Coordinate **, ViewProperties...> evaluation_points;

    //! View indicating the local id of the objects in which to evaluate the
    //! points.
    Kokkos::View<LocalOrdinal *, ViewProperties...> object_ids;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_EVALUATIONSET_HPP

//---------------------------------------------------------------------------//
// end DTK_EvaluationSet.hpp
//---------------------------------------------------------------------------//
