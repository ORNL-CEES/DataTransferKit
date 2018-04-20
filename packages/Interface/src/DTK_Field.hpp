/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file DTK_Field.hpp
 * \brief Node list.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELD_HPP
#define DTK_FIELD_HPP

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class FieldData.
 *
 * \brief Trivially-copyable field.
 *
 * \tparam Scalar The scalar type of the field degrees-of-freedom.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 */
template <class SC, class... ViewProperties>
class Field
{
  public:
    //! Scalar field type.
    using Scalar = SC;

    //! View traits.
    using ViewTraits = typename Kokkos::ViewTraits<Scalar, ViewProperties...>;

    //! The field degrees of freedom. The dof values should directly correlate
    //! to the global_dof_ids view in the dof id map. This view is rank-2 and
    //! should be dimensioned (degree of freedom, field dimension). The length
    //! of the first dimension in this view should be the same as the
    //! global_dof_ids view in the dof id map. The second dimension indicates
    //! an arbitrary field dimension. This allows for scalars, vectors, and
    //! tensors to be assigned as degrees of freedom and transferred
    //! simultaneously.
    Kokkos::View<Scalar **, ViewProperties...> dofs;
};

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_FIELD_HPP

//---------------------------------------------------------------------------//
// end DTK_Field.hpp
//---------------------------------------------------------------------------//
