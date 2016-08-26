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
 * \brief DTK_ReferenceHexMesh.hpp
 * \author Stuart R. Slattery
 * \brief Reference implementation of a parallel hex mesh for unit testing.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_REFERENCEHEXMESH_HPP
#define DTK_REFERENCEHEXMESH_HPP

#include <DTK_FunctionSpace.hpp>
#include <DTK_Field.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
namespace UnitTest
{
//---------------------------------------------------------------------------//
/*!
  \class ReferenceHexMesh
  \brief Reference parallel hex mesh implementation for unit testing.

  Provides a basic implementation of a hex mesh for unit testing. Mesh will be
  partitioned in z.
*/
//---------------------------------------------------------------------------//
class ReferenceHexMesh
{
  public:

    /*!
     * \brief Num cells constructor.
     */
    ReferenceHexMesh( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                      double x_min,
                      double x_max,
                      const int x_num_cells,
                      double y_min,
                      double y_max,
                      const int y_num_cells,
                      double z_min,
                      double z_max,
                      const int z_num_cells );

    /*!
     * \brief Edge array constructor.
     */
    ReferenceHexMesh( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                      const Teuchos::Array<double>& x_edges,
                      const Teuchos::Array<double>& y_edges,
                      const Teuchos::Array<double>& z_edges);

    /*!
     * \brief Get the function space.
     */
    Teuchos::RCP<DataTransferKit::FunctionSpace> functionSpace() const;

    /*!
     * \brief Create a field over the locally-owned nodes of the mesh.
     */
    Teuchos::RCP<DataTransferKit::Field>
    nodalField( const int field_dim ) const;

    /*!
     * \brief Create a field over the locally-owned + ghosted nodes of the
     * mesh.
     */
    Teuchos::RCP<DataTransferKit::Field>
    ghostedNodalField( const int field_dim ) const;

  private:

    // Build the mesh.
    void buildMesh( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                    const Teuchos::Array<double>& x_edges,
                    const Teuchos::Array<double>& y_edges,
                    const Teuchos::Array<double>& z_edges );

    // Build an edge array.
    Teuchos::Array<double> buildEdgeArray( 
        const double min, const double max, const int num_cells ) const;
    
    // Create a ghosted or locally-owned field.
    Teuchos::RCP<DataTransferKit::Field>
    createNodalField( const int field_dim, const bool is_ghosted ) const;

  private:

    // Function space.
    Teuchos::RCP<DataTransferKit::FunctionSpace> d_function_space;
};

//---------------------------------------------------------------------------//
} // end namespace UnitTest
} // end namespace DataTransferKit

#endif // end DTK_REFERENCEHEXMESH_HPP

//---------------------------------------------------------------------------//
// end DTK_ReferenceHexMesh.hpp
//---------------------------------------------------------------------------//
