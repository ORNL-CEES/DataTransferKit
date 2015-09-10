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
 * \brief DTK_STKMeshManager.hpp
 * \author Stuart R. Slattery
 * \brief High-level manager for STK mesh.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHMANAGER_HPP
#define DTK_STKMESHMANAGER_HPP

#include <string>
#include <unordered_map>

#include "DTK_ClientManager.hpp"
#include "DTK_Types.hpp"
#include "DTK_Field.hpp"
#include "DTK_FunctionSpace.hpp"
#include "DTK_FieldMultiVector.hpp"
#include "DTK_DBC.hpp"
#include "DTK_STKMeshField.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshManager
  \brief High-level manager for STK mesh.

  This manager provides a high-level class for automated construction of DTK
  interface objects. A user is not required to use this class but rather could
  use it to reduce code for certain implementations.
*/
//---------------------------------------------------------------------------//
class STKMeshManager : public ClientManager
{
  public:

    /*!
     * \brief Basis type enum.
     */
    enum BasisType
    {
	BASIS_TYPE_GRADIENT
    };
    
  public:

    /*!
     * \brief Default constructor.
     *
     * \param bulk_data STK mesh bulk data.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped. 
     *
     * \param basis_type The type of basis function space to use.
     */
    STKMeshManager( const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
		    const BasisType basis_type = BASIS_TYPE_GRADIENT );

    /*!
     * \brief Part name constructor.
     *
     * \param bulk_data STK mesh bulk data.
     *
     * \param part_names The names of the parts in the STK mesh that will be
     * mapped.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped. 
     *
     * \param basis_type The type of basis function space to use.
     */
    STKMeshManager( const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
		    const Teuchos::Array<std::string>& part_names,
		    const BasisType basis_type = BASIS_TYPE_GRADIENT );

    /*!
     * \brief Part vector constructor.
     *
     * \param bulk_data STK mesh bulk data.
     *
     * \param part_vector The parts in the STK mesh that will be
     * mapped.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped. 
     *
     * \param basis_type The type of basis function space to use.
     */
    STKMeshManager( const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
		    const stk::mesh::PartVector& parts,
		    const BasisType basis_type = BASIS_TYPE_GRADIENT );

    /*!
     * \brief Selector constructor.
     *
     * \param bulk_data STK mesh bulk data.
     *
     * \param selector Selector fo the parts in the STK mesh that will be
     * mapped.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped. 
     *
     * \param basis_type The type of basis function space to use.
     */
    STKMeshManager( const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
		    const stk::mesh::Selector& selector,
		    const BasisType basis_type = BASIS_TYPE_GRADIENT );

    /*!
     * \brief Register a tag with the manager that will be available for
     * solution transfer.
     */
    template<class FieldType>
    void registerField( const Teuchos::Ptr<FieldType>& field,
			const int field_dim );
    
    /*!
     * \brief Get the function space over which the mesh and its fields are
     * defined. 
     */
    Teuchos::RCP<FunctionSpace> functionSpace() const;

    /*!
     * \brief Given a field and dimension, build a vector over that field.
     */
    template<class FieldType>
    Teuchos::RCP<FieldMultiVector>
    createFieldMultiVector( const Teuchos::Ptr<FieldType>& field,
			    const int field_dim );

    //@{
    //! ClientManager interface implementation.
    /*!
     * \brief Get the entity set over which the fields are defined.
     */
    Teuchos::RCP<EntitySet> entitySet() const override;

    /*!
     * \brief Get the local map for entities supporting the function.
     */
    Teuchos::RCP<EntityLocalMap> localMap() const override;

    /*!
     * \brief Get the shape function for entities supporting the function.
     */
    Teuchos::RCP<EntityShapeFunction> shapeFunction() const override;

    /*!
     * \brief Get the integration rule for entities supporting the function.
     */
    Teuchos::RCP<EntityIntegrationRule> integrationRule() const override;

    /*!
     * \brief Get the selector function.
     */
    PredicateFunction selectFunction() const override;

    /*!
     * \brief Get the field for the given string key.
     */
    Teuchos::RCP<Field>
    field( const std::string& field_name ) const override;
    //@}

  private:

    // Create the function space.
    void createFunctionSpace( 
	const BasisType basis_type,
	const PredicateFunction& select_function );

  private:

    // Bulk data.
    Teuchos::RCP<stk::mesh::BulkData> d_bulk_data;
    
    // The function space over which the mesh and its fields are defined.
    Teuchos::RCP<FunctionSpace> d_function_space;

    // Field name to local id indexer.
    std::unordered_map<std::string,int> d_field_indexer;
    
    // Registered fields.
    Teuchos::Array<Teuchos::RCP<Field> > d_fields;
};

//---------------------------------------------------------------------------//
// Template functions.
//---------------------------------------------------------------------------//
template<class FieldType>
Teuchos::RCP<FieldMultiVector>
STKMeshManager::createFieldMultiVector( const Teuchos::Ptr<FieldType>& field,
					const int field_dim )
{
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    DTK_REQUIRE( Teuchos::nonnull(d_function_space) );
    
    Teuchos::RCP<Field> stk_field = Teuchos::rcp(
	new STKMeshField<double,FieldType>(d_bulk_data,field,field_dim) );
    return Teuchos::rcp(
	new FieldMultiVector(stk_field,d_function_space->entitySet()) );
}

//---------------------------------------------------------------------------//
template<class FieldType>
void STKMeshManager::registerField( const Teuchos::Ptr<FieldType>& field,
				    const int field_dim )
{
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
    DTK_REQUIRE( Teuchos::nonnull(d_function_space) );

    d_field_indexer.emplace( field->name(), d_fields.size() );
    
    d_fields.push_back(
	Teuchos::rcp(
	    new STKMeshField<double,FieldType>(d_bulk_data,field,field_dim) )
	);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshManager.hpp
//---------------------------------------------------------------------------//
