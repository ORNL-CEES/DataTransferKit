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

    //! View tratis.
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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_FIELD_HPP

//---------------------------------------------------------------------------//
// end DTK_Field.hpp
//---------------------------------------------------------------------------//
