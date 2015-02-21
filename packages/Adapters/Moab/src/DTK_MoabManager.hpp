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
 * \brief DTK_MoabManager.hpp
 * \author Stuart R. Slattery
 * \brief High-level manager for Moab mesh.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABMANAGER_HPP
#define DTK_MOABMANAGER_HPP

#include <string>

#include "DTK_Types.hpp"
#include "DTK_FunctionSpace.hpp"
#include "DTK_EntitySelector.hpp"

#include <Teuchos_RCP.hpp>

#include <MBParallelComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabManager
  \brief High-level manager for Moab mesh.

  This manager provides a high-level class for automated construction of DTK
  interface objects. A user is not required to use this class but rather could
  use it to reduce code for certain implementations.
*/
//---------------------------------------------------------------------------//
class MoabManager
{
  public:

    /*!
     * \brief Default constructor.
     *
     * \param moab_mesh Moab mesh.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped. 
     */
    MoabManager( const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
		 const EntityType entity_type );

    /*!
     * \brief Mesh set constructor.
     *
     * \param moab_mesh Moab mesh.
     *
     * \param mesh_set The set over which the manager will be constructed.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped. 
     */
    MoabManager( const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
		 const moab::EntityHandle& mesh_set,
		 const EntityType entity_type );

    /*!
     * \brief Destructor.
     */
    ~MoabManager();

    /*!
     * \brief Get the function space over which the mesh and its fields are
     * defined. 
     */
    Teuchos::RCP<FunctionSpace> functionSpace() const;

  private:

    // The function space over which the mesh and its fields are defined.
    Teuchos::RCP<FunctionSpace> d_function_space;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabManager.hpp
//---------------------------------------------------------------------------//
