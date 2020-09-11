/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file DTK_DOFMap.hpp
 * \brief Degree-of-freedom map.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DOFMAP_HPP
#define DTK_DOFMAP_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class DOFMap.
 *
 * \brief Trivially-copyable degree-of-freedom map.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 */
template <class... ViewProperties>
class DOFMap
{
  public:
    //! View Traits.
    using ViewTraits = typename Kokkos::ViewTraits<int, ViewProperties...>;

    //! Globally-unique ids for dofs represented on this process. These may or
    //! may not be locally owned but every dof for every object defined on
    //! this process must be available in this list. This list is of rank-1
    //! and of length equal to the number of degrees of freedom on the local
    //! MPI rank. Dimensions: (dof)
    Kokkos::View<GlobalOrdinal *, ViewProperties...> global_dof_ids;

    //! For every object of the given type in the object list give the local
    //! dof ids for that object. The local dof ids correspond to the index of
    //! the entry in the global dof id view. This view can be either rank-1 or
    //! rank-2
    //!
    //! If this view is defined as rank-1 it represents unstructured rank-2
    //! data. It should be sized as (total sum of the number of dofs defined on
    //! each object) or the total sum of the entries in the dof_per_object
    //! view. Consider the \f$n^th\f$ dof of object \f$i\f$ to be \f$d^i_n\f$
    //! which is
    //! equal to the local index of the corresponding node in the nodes
    //! view. Two objects, the first with 5 dofs and the second with 4 would
    //! then be defined via this view as: \f$(d^1_1, d^1_2, d^1_3, d^1_4, d^1_5,
    //! d^2_1, d^2_2, d^2_3, d^2_4 )\f$ with the dofs_per_object view reading
    //! \f$(5, 4)\f$. Dimensions: (object * dof).
    //!
    //! If this view is rank-2, then it represents a fixed number of dofs per
    //! object. The ordering of the dofs for each object should be the same as
    //! the case of rank-1 input. Dimensions: (object, dof).
    Kokkos::DynRankView<LocalOrdinal, ViewProperties...> object_dof_ids;

    //! The number of degrees of freedom on each object. This view is only
    //! necessary if the object_dof_ids array is rank-1. This view is rank-1
    //! and of length of the number of objects of the given type in the
    //! list. Dimensions: (object).
    Kokkos::View<unsigned *, ViewProperties...> dofs_per_object;
};

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_DOFMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_DOFMap.hpp
//---------------------------------------------------------------------------//
