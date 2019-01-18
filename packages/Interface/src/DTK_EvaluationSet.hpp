/****************************************************************************
 * Copyright (c) 2012-2019 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file DTK_EvaluationSet.hpp
 * \brief Evaluation point set.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EVALUATIONSET_HPP
#define DTK_EVALUATIONSET_HPP

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

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_EVALUATIONSET_HPP

//---------------------------------------------------------------------------//
// end DTK_EvaluationSet.hpp
//---------------------------------------------------------------------------//
