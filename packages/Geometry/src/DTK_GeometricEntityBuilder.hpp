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
 * \brief DTK_GeometricEntityBuilder.hpp
 * \author Stuart R. Slattery
 * \brief Builder for GeometricEntity classes.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRICENTITYBUILDER_HPP
#define DTK_GEOMETRICENTITYBUILDER_HPP

#include "DTK_GeometricEntity.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_AbstractFactory.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class GeometricEntityBuilder
  \brief Builder for constructing derived classes of the GeometricEntity
  class.
*/
//---------------------------------------------------------------------------//
class GeometricEntityBuilder
{
  public:

    /*!
     * \brief Constructor.
     */
    GeometricEntityBuilder();

    /*!
     * \brief Destructor.
     */
    ~GeometricEntityBuilder();

    /*!
     * \brief Set a new GeometricEntity factory with the builder.
     * \param entity_factory A factory that can create a GeometricEntity from
     * a derived class.
     * \param name A name for the factory. This should be indicative of the
     * derived class and equivalent to the entityType() field from the
     * implementation.
     */
    void setGeometricEntityFactory(
	const Teuchos::RCP<const Teuchos::AbstractFactory<GeometricEntity> >&
	entity_factory,
	const std::string& name );

    /*!
     * \brief Create a new GeometricEntity with the given name.
     * \param name Create a GeometricEntity using the factory with this name.
     * \return The created GeometricEntity.
     */
    Teuchos::RCP<GeometricEntity>
    createGeometricEntity( const std::string& name ) const;

  private:

    // Factories.
    Teuchos::Array<Teuchos::RCP<const Teuchos::AbstractFactory<GeometricEntity> > >
    d_factories;

    // Names.
    Teuchos::Array<std::string> d_names;

    // Validator.
    Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > d_validator;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_GEOMETRICENTITYBUILDER_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometricEntityBuilder.hpp
//---------------------------------------------------------------------------//
