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
 * \brief DTK_AbstractObjectRegistry.hpp
 * \author Stuart R. Slattery
 * \brief Derived class registration with base class.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ABSTRACTOBJECTREGISTRY_HPP
#define DTK_ABSTRACTOBJECTREGISTRY_HPP

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class UndefinedAbstractObjectRegistrationPolicy
  \brief Complie time error indicator for policy implementations.
*/
//---------------------------------------------------------------------------//
template<typename Derived>
struct UndefinedAbstractObjectRegistrationPolicy 
{
    static inline Derived notDefined() 
    {
	return Derived::this_type_is_missing_a_specialization();
    }
};

//---------------------------------------------------------------------------//
/*!
  \class AbstractObjectRegistrationPolicy
  \brief Policy definition for objects that can be built through an
  abstract builder.
*/
//---------------------------------------------------------------------------//
template<class Derived>
class AbstractObjectRegistrationPolicy
{
  public:

    //! Base class type.
    typedef Derived object_type;

    /*!
     * \brief Register a derived class with a base class.
     */
    static void registerDerivedClassWithBaseClass()
    {
	UndefinedAbstractObjectRegistrationPolicy<Derived>::notDefined();
	return Teuchos::null;
    }
};

//---------------------------------------------------------------------------//
// Specialization for void.
template<>
class AbstractObjectRegistrationPolicy<void>
{
  public:
    typedef void object_type;
    static void registerDerivedClassWithBaseClass()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
/*!
  \class AbstractObjectRegistry
  \brief Compile-time registy from derived objects inheriting from a base
  class.

  This class lets the user indicate which derived types will be used with a
  given base class. Derived classes must be registered with base classes
  inheriting from abstract compile time interfaces before use.
*/
//---------------------------------------------------------------------------//
template<class Base, 
	 class Derived1 = void, 
	 class Derived2 = void, 
	 class Derived3 = void,
	 class Derived4 = void, 
	 class Derived5 = void, 
	 class Derived6 = void,
	 class Derived7 = void, 
	 class Derived8 = void, 
	 class Derived9 = void>
class AbstractObjectRegistry
{
  public:

    //! Typedefs.
    typedef Base     base_type;
    typedef Derived1 derived_type_1;
    typedef Derived1 derived_type_2;
    typedef Derived1 derived_type_3;
    typedef Derived1 derived_type_4;
    typedef Derived1 derived_type_5;
    typedef Derived1 derived_type_6;
    typedef Derived1 derived_type_7;
    typedef Derived1 derived_type_8;
    typedef Derived1 derived_type_9;

    /*!
     * \brief Constructor.
     */
    AbstractObjectRegistry()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~AbstractObjectRegistry()
    { /* ... */ }

    /*!
     * \brief Set an abstract factory for a AbstractObjectRegistry subclass.
     * \param factory A factory for a AbstractObjectRegistry subclass.
     */
    static void registerDerivedClasses();
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_AbstractObjectRegistry_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_ABSTRACTOBJECTREGISTRY_HPP

//---------------------------------------------------------------------------//
// end DTK_AbstractObjectRegistry.hpp
//---------------------------------------------------------------------------//
