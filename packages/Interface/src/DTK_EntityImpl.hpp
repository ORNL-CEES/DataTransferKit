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
 * \brief DTK_EntityImpl.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYIMPL_HPP
#define DTK_ENTITYIMPL_HPP

#include <string>

#include "DTK_Types.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityImpl
  \brief Geometric entity implementation definition.
*/
//---------------------------------------------------------------------------//
class EntityImpl
{
  public:

    /*!
     * \brief Constructor.
     */
    EntityImpl();

    /*!
     * \brief Destructor.
     */
    virtual ~EntityImpl();

    //@{
    //! Identification functions.
    /*!
     * \brief Return a string indicating the derived entity type.
     * \return A string indicating the type of derived entity implementing the
     * interface.
     */
    virtual std::string name() const;

    /*!
     * \brief Get the entity type.
     * \return The entity type.
     */
    virtual EntityType entityType() const;

    /*!
     * \brief Get the unique global identifier for the entity.
     * \return A unique global identifier for the entity.
     */
    virtual EntityId id() const;
    
    /*!
     * \brief Get the parallel rank that owns the entity.
     * \return The parallel rank that owns the entity.
     */
    virtual int ownerRank() const;
    //@}

    //@{
    //! Geometric functions.
    /*!
     * \brief Return the physical dimension of the entity.
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    virtual int physicalDimension() const;

    /*!
     * \brief Return the axis-aligned bounding box around the entity.
     * \param bounding_box A Cartesian box that bounds the entity.
     */
    virtual void boundingBox( Teuchos::Tuple<double,6>& bounds ) const;
    //@}

    //@{
    //! AbstractBuildableObjectPolicy interface.
    /*!
     * \brief Return a string indicating the derived object type.
     * \return A string indicating the type of derived object implementing the
     * interface.
     */
    std::string objectType() const;
    //@}

    //@{
    //! AbstractSerializableObjectPolicy interface.
    /*
     * \brief Serialize the entity into a buffer.
     * \param buffer A view into a data buffer of size byteSize(). Write the
     * serialized entity into this view.
     */
    virtual void serialize( const Teuchos::ArrayView<char>& buffer ) const;

    /*!
     * \brief Deserialize an entity from a buffer.
     * \param buffer A view into a data buffer of size byteSize(). Deserialize
     * the entity from this view.
     */
    virtual void deserialize( const Teuchos::ArrayView<const char>& buffer );
    //@}
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYIMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityImpl.hpp
//---------------------------------------------------------------------------//
