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
 * \brief DTK_STKMeshField.hpp
 * \author Stuart R. Slattery
 * \brief STK field data access.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHFIELD_HPP
#define DTK_STKMESHFIELD_HPP

#include <vector>
#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_Field.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ArrayView.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshField
  \brief Field data access for STK mesh.
*/
//---------------------------------------------------------------------------//
template<class Scalar,class FieldType>
class STKMeshField : public Field
{
  public:

    /*!
     * \brief Constructor.
     */
    STKMeshField( const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
                  const Teuchos::Ptr<FieldType>& field,
                  const int field_dim );

    /*!
     * \brief Get the dimension of the field.
     */
    int dimension() const override;

    /*!
     * \brief Get the locally-owned entity support location ids of the field.
     */
    Teuchos::ArrayView<const SupportId> getLocalSupportIds() const override;

    /*!
     * \brief Given a local support id and a dimension, read data from the
     * application field.
     */
    double readFieldData( const SupportId support_id,
                          const int dimension ) const override;

    /*!
     * \brief Given a local support id, dimension, and field value, write data
     * into the application field.
     */
    void writeFieldData( const SupportId support_id,
                         const int dimension,
                         const double data ) override;

    /*!
     * \brief Finalize a field after writing into it.
     */
    void finalizeAfterWrite() override;

  private:

    // The mesh over which the field is defined.
    Teuchos::RCP<stk::mesh::BulkData> d_bulk_data;

    // The field containing the vector data.
    Teuchos::Ptr<FieldType> d_field;

    // The dimension of the field.
    int d_field_dim;

    // The enitities over which the field is defined.
    std::vector<stk::mesh::Entity> d_field_entities;

    // The support ids of the entities over which the field is constructed.
    Teuchos::Array<SupportId> d_support_ids;

    // Support id to local id map.
    std::unordered_map<SupportId,int> d_id_map;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_STKMeshField_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHFIELD_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshField.hpp
//---------------------------------------------------------------------------//
