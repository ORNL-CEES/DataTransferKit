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
 * \brief DTK_FieldGroup.hpp
 * \author Stuart R. Slattery
 * \brief Field group.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDGROUP_HPP
#define DTK_FIELDGROUP_HPP

#include <string>

#include "DTK_EntitySet.hpp"
#include "DTK_EntityFunctionSpace.hpp"
#include "DTK_DOFAccessor.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class FieldGroup
  \brief Group for aggregating entity sets, fields, and their functional
  support.

  The FieldGroup binds a particular discretization of a field
  (EntityFunctionSpace) to an entity set. Multiple fields can be supported by
  this discretization and access to the degrees of freedom for each of these
  fields is managed by the DOFAccessor. The FieldGroup lets clients register
  multiple degree of freedom sets with a particular discretization on an
  entity set.
*/
//---------------------------------------------------------------------------//
class FieldGroup
{
  public:

    /*!
     * \brief Constructor.
     * \param entity_set The entity set over which the fields are defined.
     * \param entity_space Accessor to the functional support of the fields.
     */
    FieldGroup( const Teuchos::RCP<EntitySet>& entity_set,
		const Teuchos::RCP<EntityFunctionSpace>& entity_space );

    /*!
     * \brief Destructor.
     */
    ~FieldGroup();

    /*!
     * \brief Add a DOF accessor for a given field to the group.
     * \param dof_accessor Accessor to the degrees of freedom for a field
     * supported by this entity set and functional support.
     * \param name The name of the field the degrees of freedom represent.
     */
    void registerDOFAccessor( const Teuchos::RCP<DOFAccessor>& dof_accessor,
			      const std::string& name );

    /*!
     * \brief Get the entity set over which the fields are defined.
     */
    Teuchos::RCP<EntitySet> entitySet() const;

    /*!
     * \brief Get the functional support accessor for the fields.
     */
    Teuchos::RCP<EntityFunctionSpace> entitySpace() const;

    /*!
     * \brief Get a DOF accessor with a given name.
     * \param The name of the field the degrees of freedom represent.
     * \return Accessor to the degrees of freedom for a field supported by
     * this entity set and functional support.
     */
    Teuchos::RCP<DOFAccessor> getDOFAccessor( const std::string& name );

  private:

    // The entity set over which the field is constructed.
    Teuchos::RCP<EntitySet> d_entity_set;

    // The field support space for the entities in the set.
    Teuchos::RCP<EntityFunctionSpace> d_entity_space;

    // DOF Accessors.
    Teuchos::Array<Teuchos::RCP<DOFAccessor> > d_dof_accessor;

    // DOF Accessor names.
    Teuchos::Array<std::string> d_accessor_names;

    // DOF Accessor validator.
    Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    d_accessor_validator;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_FIELDGROUP_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldGroup.hpp
//---------------------------------------------------------------------------//
