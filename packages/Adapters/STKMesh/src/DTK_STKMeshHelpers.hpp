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
 * \brief DTK_STKMeshHelpers.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHHELPERS_HPP
#define DTK_STKMESHHELPERS_HPP

#include <string>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Intrepid_FieldContainer.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshHelpers
  \brief A stateless class of helpers for STK mesh.
*/
//---------------------------------------------------------------------------//
class STKMeshHelpers
{
  public:

    /*!
     * \brief Constructor.
     */
    STKMeshHelpers() { /* ... */ } 

    /*!
     * \brief Destructor.
     */
    ~STKMeshHelpers() { /* ... */ }

    /*!
     * \brief Given a DTK entity, extract the STK entity.
     */
    static const stk::mesh::Entity& extractEntity( const Entity dtk_entity );

    /*!
     * \brief Given a DTK EntityType, get the STK entity rank.
     */
    static stk::mesh::EntityRank getRankFromType( 
	const EntityType dtk_type, const int space_dim );

    /*!
     * \brief Given a STK entity rank, get the DTK entity type.
     */
    static EntityType getTypeFromRank( 
	const stk::mesh::EntityRank stk_rank, const int space_dim );

    /*!
     * \brief Given a DTK entity, return the corresponding STK entity key.
     */
    static stk::mesh::EntityKey getKeyFromEntity( const Entity dtk_entity );

    /*!
     * \brief Given a STK entity, return the coordinates of its nodes in a
     * field container ordered by canonical node order (C,N,D).
     */
    static Intrepid::FieldContainer<double> 
    getEntityNodeCoordinates( 
	const Teuchos::Array<stk::mesh::Entity>& stk_entities,
	const stk::mesh::BulkData& bulk_data );

  private:

    /*!
     * \brief Given a STK entity, return the coordinates of its nodes in a
     * field container ordered by canonical node order (C,N,D). Templated
     * extraction layer.
     */
    template<class FieldType>
    static Intrepid::FieldContainer<double> 
    extractEntityNodeCoordinates( 
	const Teuchos::Array<stk::mesh::Entity>& stk_entities,
	const stk::mesh::BulkData& bulk_data,
	const int space_dim );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_STKMeshHelpers_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHHELPERS_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshHelpers.hpp
//---------------------------------------------------------------------------//
