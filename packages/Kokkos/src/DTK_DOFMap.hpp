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
 * \file DTK_DOFMap.hpp
 * \brief Degree-of-freedom map.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DOFMAP_HPP
#define DTK_DOFMAP_HPP

#include "DTK_ConfigDefs.hpp"

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
    //! may not be equally owned but every dof for every object defined on this
    //! process (ghosted or not) must be available in this list. This list is
    //! of rank-1 and of length equal to the number of degrees of freedom on
    //! the local MPI rank. Dimensions: (dof)
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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_DOFMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_DOFMap.hpp
//---------------------------------------------------------------------------//
