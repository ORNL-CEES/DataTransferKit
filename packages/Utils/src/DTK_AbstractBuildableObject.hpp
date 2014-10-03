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
 * \brief DTK_AbstractBuildableObject.hpp
 * \author Stuart R. Slattery
 * \brief Abstract buildable object interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ABSTRACTBUILDABLEOBJECT_HPP
#define DTK_ABSTRACTBUILDABLEOBJECT_HPP

#include <string>

#include "DTK_AbstractBuilder.hpp"
#include "DTK_AbstractSerializableObjectPolicy.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class UndefinedAbstractBuildableObjectPolicy
  \brief Complie time error indicator for policy implementations.
*/
//---------------------------------------------------------------------------//
template<typename T>
struct UndefinedAbstractBuildableObjectPolicy 
{
    static inline T notDefined() 
    {
	return T::this_type_is_missing_a_specialization();
    }
};

//---------------------------------------------------------------------------//
/*!
  \class AbstractBuildableObjectPolicy
  \brief Policy definition for objects that can be built through an
  abstract builder.

*/
//---------------------------------------------------------------------------//
template<class T>
class AbstractBuildableObjectPolicy
{
  public:

    //! Base class type.
    typedef T object_type;

    //@{
    //! Identification functions.
    /*!
     * \brief Return a string indicating the derived object type.
     * \return A string indicating the type of derived object implementing the
     * interface. This string will drive object construction with the builder.
     */
    static std::string objectType( const Teuchos::RCP<T>& object )
    {
	UndefinedAbstractBuildableObjectPolicy<T>::notDefined();
	return std::string("Not implemented");
    }
    //@}

    //@{
    //! Polymorphic construction functions.
    /*!
     * \brief Static function for getting the builder for the base class.
     */
    static Teuchos::RCP<AbstractBuilder<T> > getBuilder()
    {
	UndefinedAbstractSerializableObjectPolicy<T>::notDefined();
	return Teuchos::null;
    }
    //@}
};

//---------------------------------------------------------------------------//
/*!
  \class AbstractBuildableObject
  \brief Interface definition for objects that can be built through an
  abstract builder.

  This class provides the ability to attach a static builder to the derived
  type. Users can then register factories for subclasses of the derived type
  and access the builder when needed.
*/
//---------------------------------------------------------------------------//
template<class Object>
class AbstractBuildableObject
{
  public:

    //! Typedefs.
    typedef AbstractBuildableObjectPolicy<Object> ABOP;

    /*!
     * \brief Constructor.
     */
    AbstractBuildableObject();

    /*!
     * \brief Destructor.
     */
    virtual ~AbstractBuildableObject();

    /*!
     * \brief Set an abstract factor for a AbstractBuildableObject subclass.
     * \param builder A factory for a AbstractBuildableObject subclass.
     */
    static void setDerivedClassFactory(
	const Teuchos::RCP<const Teuchos::AbstractFactory<Object> >& factory );

    /*!
     * \brief Get the abstract builder for AbstractBuildableObject subclasses.
     * \return The builder for AbstractBuildableObject subclasses.
     */
    static Teuchos::RCP<AbstractBuilder<Object> > getBuilder();
    //@}

  private:
    
    // Abstract builder for geometric entity subclasses.
    static Teuchos::RCP<AbstractBuilder<Object> > b_builder;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_AbstractBuildableObject_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_ABSTRACTBUILDABLEOBJECT_HPP

//---------------------------------------------------------------------------//
// end DTK_AbstractBuildableObject.hpp
//---------------------------------------------------------------------------//
