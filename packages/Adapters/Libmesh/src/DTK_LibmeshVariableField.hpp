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
 * \brief DTK_LibmeshVariableField.hpp
 * \author Stuart R. Slattery
 * \brief Libmesh system vector data access.
 */
//---------------------------------------------------------------------------//

#ifndef LIBMESHDTKADAPTERS_LIBMESHVARIABLEFIELD_HPP
#define LIBMESHDTKADAPTERS_LIBMESHVARIABLEFIELD_HPP

#include <string>
#include <unordered_map>

#include <DTK_Field.hpp>
#include <DTK_Types.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <libmesh/mesh_base.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class LibmeshVariableField
  \brief DTK field implementation for libmesh variables.
*/
//---------------------------------------------------------------------------//
class LibmeshVariableField : public Field
{
  public:
    /*!
     * \brief Constructor.
     * \param libmesh_mesh The mesh.
     * \param libmesh_system The system containing the variables.
     * \param variable_name The name of the variable for which we will
     * create the vector. The vector will be defined over all active
     * subdomains of this variable.
     */
    LibmeshVariableField( const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
                          const Teuchos::RCP<libMesh::System> &libmesh_system,
                          const std::string &variable_name );

    /*!
     * \brief Get the dimension of the field.
     */
    int dimension() const;

    /*!
     * \brief Get the locally-owned entity support ids of the field.
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
    void writeFieldData( const SupportId support_id, const int dimension,
                         const double data ) override;

    /*!
     * \brief Finalize writing of field data to a field. This lets some
     * clients do a write post-process (e.g. update ghost values).
     */
    void finalizeAfterWrite() override;

  private:
    // Libmesh mesh.
    Teuchos::RCP<libMesh::MeshBase> d_libmesh_mesh;

    // Libmesh system.
    Teuchos::RCP<libMesh::System> d_libmesh_system;

    // System id.
    int d_system_id;

    // Variable id.
    int d_variable_id;

    // The support ids of the entities over which the field is constructed.
    Teuchos::Array<SupportId> d_support_ids;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end LIBMESHDTKADAPTERS_LIBMESHVARIABLEFIELD_HPP

//---------------------------------------------------------------------------//
// end DTK_LibmeshVariableField.hpp
//---------------------------------------------------------------------------//
