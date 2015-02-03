//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_MoabMeshSetIndexer.hpp
 * \author Stuart R. Slattery
 * \brief Moab mesh set indexer.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABMESHSETINDEXER_HPP
#define DTK_MOABMESHSETINDEXER_HPP

#include <unordered_map>

#include <Teuchos_RCP.hpp>

#include <MBParallelComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabMeshSetIndexer
  \brief Moab mesh set indexer.
*/
//---------------------------------------------------------------------------//
class MoabMeshSetIndexer
{
  public:

    /*!
     * \brief Constructor.
     * \param mesh The moab interface.
     */
    MoabMeshSetIndexer( const Teuchos::RCP<moab::ParallelComm>& moab_mesh );

    /*!
     * \brief Destructor.
     */
     ~MoabMeshSetIndexer();

    /*!
     * \brief Given an entity set handle, get the integer index in the mesh.
     */
    int getIndexFromMeshSet( const moab::EntityHandle mesh_set ) const;

    /*!
     * \brief Given an integer index, get the entity set handle.
     */
    moab::EntityHandle getMeshSetFromIndex( const int index ) const;

  private:

    // Moab mesh.
    Teuchos::RCP<moab::ParallelComm> d_moab_mesh;

    // Handle-to-index map.
    std::unordered_map<moab::EntityHandle,int> d_handle_to_index_map;

    // Index-to-handle map.
    std::unordered_map<int,moab::EntityHandle> d_index_to_handle_map;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABMESHSETINDEXER_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabMeshSetIndexer.hpp
//---------------------------------------------------------------------------//
