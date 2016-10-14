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
 * \brief DTK_EntityCenteredField.hpp
 * \author Stuart R. Slattery
 * \brief Entity-centered field.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYCENTEREDFIELD_HPP
#define DTK_ENTITYCENTEREDFIELD_HPP

#include <unordered_map>

#include "DTK_Entity.hpp"
#include "DTK_Field.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityCenteredField
  \brief Entity-centered field.

  Field implementation for basic geometry fields over entities. Input data is
  1D but ordered in 2D such that entity 'e' has for dimension 'd' a value of:

  BLOCKED: data[d][e]

  INTERLEAVED: data[e][d]

  Set the value of data layout in the constructor.
*/
//---------------------------------------------------------------------------//
class EntityCenteredField : public Field
{
  public:
    //! Blocked/Interleaved data layout enum.
    enum DataLayout
    {
        BLOCKED,
        INTERLEAVED
    };

    /*!
     * \brief Entity constructor.
     */
    EntityCenteredField( const Teuchos::ArrayView<Entity> &entities,
                         const int field_dim,
                         const Teuchos::ArrayRCP<double> &dof_data,
                         const DataLayout layout );

    /*!
     * \brief Entity id constructor.
     */
    EntityCenteredField( const Teuchos::ArrayView<const EntityId> &entity_ids,
                         const int field_dim,
                         const Teuchos::ArrayRCP<double> &dof_data,
                         const DataLayout layout );

    /*!
     * \brief Get the dimension of the field.
     */
    int dimension() const;

    /*!
     * \brief Get the locally-owned support location ids of the field.
     */
    Teuchos::ArrayView<const SupportId> getLocalSupportIds() const;

    /*!
     * \brief Given a local dof id and a dimension, read data from the
     * application field.
     */
    double readFieldData( const SupportId support_id,
                          const int dimension ) const;

    /*!
     * \brief Given a local dof id, dimension, and field value, write data
     * into the application field.
     */
    void writeFieldData( const SupportId support_id, const int dimension,
                         const double data );

  private:
    // The dof ids of the entities over which the field is constructed.
    Teuchos::Array<SupportId> d_support_ids;

    // The dimension of the field.
    int d_field_dim;

    // The field data.
    Teuchos::ArrayRCP<double> d_data;

    // Data layout.
    DataLayout d_layout;

    // Field leading dimension.
    int d_lda;

    // Dof id to local id map.
    std::unordered_map<SupportId, int> d_id_map;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYCENTEREDFIELD_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityCenteredField.hpp
//---------------------------------------------------------------------------//
