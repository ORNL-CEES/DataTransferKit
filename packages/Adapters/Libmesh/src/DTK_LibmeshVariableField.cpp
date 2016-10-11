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
 * \brief DTK_LibmeshVariableField.cpp
 * \author Stuart R. Slattery
 * \brief Libmesh system vector manager
 */
//---------------------------------------------------------------------------//

#include <cassert>
#include <unordered_map>

#include "DTK_DBC.hpp"
#include "DTK_LibmeshVariableField.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/node.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
LibmeshVariableField::LibmeshVariableField(
    const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::RCP<libMesh::System> &libmesh_system,
    const std::string &variable_name )
    : d_libmesh_mesh( libmesh_mesh )
    , d_libmesh_system( libmesh_system )
{
    // Get ids.
    d_system_id = d_libmesh_system->number();
    d_variable_id = d_libmesh_system->variable_number( variable_name );

    // Get the local support ids.
    libMesh::MeshBase::const_node_iterator nodes_end =
        d_libmesh_mesh->local_nodes_end();
    for ( libMesh::MeshBase::const_node_iterator node_it =
              d_libmesh_mesh->local_nodes_begin();
          node_it != nodes_end; ++node_it )
    {
        DTK_CHECK( ( *node_it )->valid_id() );
        d_support_ids.push_back( ( *node_it )->id() );
    }
}

//---------------------------------------------------------------------------//
// Get the dimension of the field.
int LibmeshVariableField::dimension() const
{
    // Libmesh variables are only scalars.
    return 1;
}

//---------------------------------------------------------------------------//
// Get the locally-owned entity support ids of the field.
Teuchos::ArrayView<const SupportId>
LibmeshVariableField::getLocalSupportIds() const
{
    return d_support_ids();
}

//---------------------------------------------------------------------------//
// Given a local support id and a dimension, read data from the application
// field.
double LibmeshVariableField::readFieldData( const SupportId support_id,
                                            const int dimension ) const
{
    DTK_REQUIRE( 0 == dimension );
    const libMesh::Node &node = d_libmesh_mesh->node( support_id );
    DTK_CHECK( 1 == node.n_comp( d_system_id, d_variable_id ) );
    libMesh::dof_id_type dof_id =
        node.dof_number( d_system_id, d_variable_id, 0 );
    return d_libmesh_system->current_local_solution->el( dof_id );
}

//---------------------------------------------------------------------------//
// Given a local support id, dimension, and field value, write data into the
// application field.
void LibmeshVariableField::writeFieldData( const SupportId support_id,
                                           const int dimension,
                                           const double data )
{
    DTK_REQUIRE( 0 == dimension );
    const libMesh::Node &node = d_libmesh_mesh->node( support_id );
    DTK_CHECK( 1 == node.n_comp( d_system_id, d_variable_id ) );
    if ( node.processor_id() == d_libmesh_system->processor_id() )
    {
        libMesh::dof_id_type dof_id =
            node.dof_number( d_system_id, d_variable_id, 0 );
        d_libmesh_system->solution->set( dof_id, data );
    }
}

//---------------------------------------------------------------------------//
// Finalize after writing.
void LibmeshVariableField::finalizeAfterWrite()
{
    d_libmesh_system->solution->close();
    d_libmesh_system->update();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshVariableField.cpp
//---------------------------------------------------------------------------//
