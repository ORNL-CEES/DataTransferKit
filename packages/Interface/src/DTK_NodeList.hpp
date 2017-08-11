/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file DTK_NodeList.hpp
 * \brief Node list.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_NODELIST_HPP
#define DTK_NODELIST_HPP

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class NodeList.
 *
 * \brief Trivially-copyable node list.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 */
template <class... ViewProperties>
class NodeList
{
  public:
    //! View Traits
    using ViewTraits = Kokkos::ViewTraits<int, ViewProperties...>;

    //! The coordinates of the nodes that are locally-owned by this MPI
    //! rank. This view is rank-2 and should be sized as (number of nodes,
    //! spatial dimension)
    Kokkos::View<Coordinate **, ViewProperties...> coordinates;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_NODELIST_HPP

//---------------------------------------------------------------------------//
// end DTK_NodeList.hpp
//---------------------------------------------------------------------------//
