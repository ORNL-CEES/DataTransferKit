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
 * \brief DTK_AbstractBuilder.cpp
 * \author Stuart R. Slattery
 * \brief Builder for Abstract classes.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ABSTRACTBUILDER_IMPL_HPP
#define DTK_ABSTRACTBUILDER_IMPL_HPP

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Base>
AbstractBuilder<Base>::AbstractBuilder()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
template<class Base>
AbstractBuilder<Base>::~AbstractBuilder()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Set a new Abstract factory with the builder.
template<class Base>
void AbstractBuilder<Base>::setDerivedClassFactory(
    const Teuchos::RCP<const Teuchos::AbstractFactory<Base> >&
    entity_factory,
    const std::string& name )
{
    d_factories.push_back( entity_factory );
    d_names.push_back( name );
}

//---------------------------------------------------------------------------//
// Get the integral key for a string key.
template<class Base>
int AbstractBuilder<Base>::getIntegralKey( const std::string& name )
{
    buildValidator();
    return d_validator->getIntegralValue( name );
}

//---------------------------------------------------------------------------//
// Create a new Abstract with the given name.
template<class Base>
Teuchos::RCP<Base>
AbstractBuilder<Base>::create( const std::string& name )
{
    buildValidator();
    int factory_index = d_validator->getIntegralValue( name );
    return d_factories[ factory_index ]->create();
}
//---------------------------------------------------------------------------//
// Create a new Abstract with the given integral key.
template<class Base>
Teuchos::RCP<Base>
AbstractBuilder<Base>::create( const int key )
{
    return d_factories[ key ]->create();
}

//---------------------------------------------------------------------------//
// Build the validator.
template<class Base>
void AbstractBuilder<Base>::buildValidator()
{
    if ( Teuchos::is_null(d_validator) )
    {
	std::string default_name( "No Default Implementation" );
	d_validator = Teuchos::rcp(
	    new Teuchos::StringToIntegralParameterEntryValidator<int>(
		d_names(), default_name) );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ABSTRACTBUILDER_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_AbstractBuilder_impl.hpp
//---------------------------------------------------------------------------//
