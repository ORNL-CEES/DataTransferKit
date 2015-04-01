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
 * \brief DTK_MoabTagField.hpp
 * \author Stuart R. Slattery
 * \brief Moab tag vector manager.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABTAGFIELD_HPP
#define DTK_MOABTAGFIELD_HPP

#include <vector>
#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_Field.hpp"

#include <MBParallelComm.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabTagField
  \brief Field implementation for moab tags.
*/
//---------------------------------------------------------------------------//
template<class Scalar>
class MoabTagField : public Field<Scalar>
{
  public:

    /*!
     * \brief Constructor.
     * \param moab_mesh The mesh containing the mesh set and tag.
     * \param mesh_set The set of mesh entities over which the vector is defined.
     * \param tag The tag containing the vector data.
     */
    MoabTagField( const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
		  const moab::EntityHandle& mesh_set,
		  const moab::Tag& tag );

    /*!
     * \brief Get the dimension of the field.
     */
    int dimension() const;

    /*!
     * \brief Get the locally-owned support location ids of the field.
     */
    Teuchos::ArrayView<const SupportId> getLocalSupportIds() const;

    /*!
     * \brief Given a local support id and a dimension, read data from the
     * application field.
     */
    Scalar readFieldData( const SupportId support_id,
			  const int dimension ) const;

    /*!
     * \brief Given a local support id, dimension, and field value, write data
     * into the application field.
     */
    void writeFieldData( const SupportId support_id,
			 const int dimension,
			 const Scalar data );

  private:

    // The mesh over which the tag is defined.
    Teuchos::RCP<moab::ParallelComm> d_moab_mesh;

    // The mesh set over which the vector is defined.
    moab::EntityHandle d_mesh_set;

    // The tag containing the vector data.
    moab::Tag d_tag;

    // The dimension of the tag.
    int d_tag_dim;

    // Entities in the mesh set.
    std::vector<moab::EntityHandle> d_entities;
    
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

#include "DTK_MoabTagField_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_MOABTAGFIELD_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabTagField.hpp
//---------------------------------------------------------------------------//
