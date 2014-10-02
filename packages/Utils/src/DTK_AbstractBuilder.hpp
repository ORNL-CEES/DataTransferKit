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
 * \brief DTK_AbstractBuilder.hpp
 * \author Stuart R. Slattery
 * \brief Builder for abstract classes.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ABSTRACTBUILDER_HPP
#define DTK_ABSTRACTBUILDER_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_AbstractFactory.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class AbstractBuilder
  \brief Builder for constructing derived classes of the base.
*/
//---------------------------------------------------------------------------//
template<class Base>
class AbstractBuilder
{
  public:

    /*!
     * \brief Constructor.
     */
    AbstractBuilder();

    /*!
     * \brief Destructor.
     */
    ~AbstractBuilder();

    /*!
     * \brief Set a new Abstract factory with the builder.
     * \param entity_factory A factory that can create a Abstract from
     * a derived class.
     * \param name A name for the factory. This should be indicative of the
     * derived class and equivalent to the entityType() field from the
     * implementation.
     */
    void setDerivedClassFactory(
	const Teuchos::RCP<const Teuchos::AbstractFactory<Base> >&
	entity_factory,
	const std::string& name );

    /*!
     * \brief Get the integral key for a string key.
     * \param The name to get the key for.
     * \return The integral key.
     */
    int getIntegralKey( const std::string& name );

    /*!
     * \brief Create a new base class with the given derived class factory
     * name.  
     * \param name Create a base class using the derived class factory with
     * this name.
     * \return The created base class object.
     */
    Teuchos::RCP<Base> create( const std::string& name );
          
    /*!
     * \brief Create a new base class with the given derived class factory
     * integral key.  
     * \param name Create a base class using the derived class factory with
     * this integral key.
     * \return The created base class object.
     */
    Teuchos::RCP<Base> create( const int key );

  private:

    // Build the validator.
    void buildValidator();

  private:

    // Factories.
    Teuchos::Array<Teuchos::RCP<const Teuchos::AbstractFactory<Base> > >
    d_factories;

    // Names.
    Teuchos::Array<std::string> d_names;

    // Validator.
    Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    d_validator;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_AbstractBuilder_impl.hpp"

#endif // end DTK_ABSTRACTBUILDER_HPP

//---------------------------------------------------------------------------//
// end DTK_AbstractBuilder.hpp
//---------------------------------------------------------------------------//
