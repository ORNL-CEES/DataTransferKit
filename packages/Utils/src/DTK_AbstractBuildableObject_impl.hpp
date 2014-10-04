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
 * \brief DTK_AbstractBuildableObject_impl.hpp
 * \author Stuart R. Slattery
 * \brief Abstract buildable object interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ABSTRACTBUILDABLEOBJECT_IMPL_HPP
#define DTK_ABSTRACTBUILDABLEOBJECT_IMPL_HPP

#include "DTK_DBC.hpp"

#include <Teuchos_AbstractFactoryStd.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Object>
AbstractBuildableObject<Object>::AbstractBuildableObject()
{ /* ... */ }

//---------------------------------------------------------------------------//
//brief Destructor.
template<class Object>
AbstractBuildableObject<Object>::~AbstractBuildableObject()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Static members.
//---------------------------------------------------------------------------//
// Add factories through the setDerivedClassFactory function.
template<class Object>
Teuchos::RCP<AbstractBuilder<Object> > 
AbstractBuildableObject<Object>::b_builder = 
    Teuchos::rcp( new AbstractBuilder<Object>() );

//---------------------------------------------------------------------------//
// Set an abstract builder for AbstractBuildableObject subclasses.
template<class Object>
template<class DerivedObject>
void AbstractBuildableObject<Object>::setDerivedClassFactory()
{ 
    Teuchos::RCP<const Teuchos::AbstractFactory<Object> > factory =
	Teuchos::abstractFactoryStd<Object,DerivedObject>();
    b_builder->setDerivedClassFactory(
	factory, ABOP::objectType(factory->create()) );
}

//---------------------------------------------------------------------------//
// Get an abstract builder for AbstractBuildableObject subclasses.
template<class Object>
Teuchos::RCP<AbstractBuilder<Object> > 
AbstractBuildableObject<Object>::getBuilder()
{ 
    return b_builder; 
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ABSTRACTBUILDABLEOBJECT_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_AbstractBuildableObject_impl.hpp
//---------------------------------------------------------------------------//
