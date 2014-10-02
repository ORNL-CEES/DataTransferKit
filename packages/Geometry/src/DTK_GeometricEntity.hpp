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
 * \brief DTK_GeometricEntity.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRICENTITY_HPP
#define DTK_GEOMETRICENTITY_HPP

#include <string>

#include "DTK_GeometryTypes.hpp"
#include "DTK_MappingStatus.hpp"
#include "DTK_AbstractBuilder.hpp"
#include "DTK_SerializableAbstractObjectPolicy.hpp"
#include "DTK_AbstractSerializer.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_SerializationTraits.hpp>

namespace DataTransferKit
{
// Forward declaration of box.
class Box;

//---------------------------------------------------------------------------//
/*!
  \class GeometricEntity
  \brief Geometric entity interface definition.

  GeometricEntity provides access to basic properites of geometric objects. A
  geometry is simply an object or collection of objects that has $n$ physical
  dimensions and a spatial domain $\Omega \in \mathbb{R}^n$ that is bounded by
  a boundary $\Gamma \in \mathbb{R}^n$. Concrete examples of geometries in 3
  dimensions include cubes, cylinders, polyhedron, or mesh elements. A
  geometry can have 1, 2, or three dimensions. To specify the general position
  in space of the geometry, each object is required to have a centroid given
  in Cartesian coordinates with (x) given for 1 dimensional geometries, (x,y)
  for two dimensions, and (x,y,z) for 3 dimensions. A measure is also
  specified for each geometry where the measure is defined as length in 1
  dimension, area in 2 dimensions, and volume for 3 dimensions. In addition to
  this data, a geometry must be able to provide a Cartesian axis-aligned
  bounding box that encapsulates the entire geometry. For geometric search
  operations to be performed, a geometry must be able to determine if a given
  point of the same dimensionality as the geometry is contained within the
  boundary of the geometry (i.e. $\hat{r} \in \Omega$).
*/
//---------------------------------------------------------------------------//
class GeometricEntity
{
  public:

    /*!
     * \brief Constructor.
     */
    GeometricEntity();

    /*!
     * \brief Destructor.
     */
    virtual ~GeometricEntity();

    //@{
    //! Identification functions.
    /*!
     * \brief Return a string indicating the derived entity type.
     * \return A string indicating the type of derived entity implementing the
     * interface.
     */
    virtual std::string entityType() const;

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
     * \brief Return the parametric dimension of the entity.
     * \return The parameteric dimension of the entity. This is the dimension
     * of coordinates needed to parametrize the space of an entity. This may
     * be different than the physical dimension.
     */
    virtual int parametricDimension() const;

    /*!
     * \brief Return the entity measure with respect to the parameteric
     * dimension (volume for a 3D entity, area for 2D, and length for 1D).
     * \return The measure of the entity.
     */
    virtual double measure() const;

    /*!
     * \brief Return the centroid of the entity.
     * \param centroid A view to an array of size physicalDimension(). Write
     * the centroid into this view.
     */
    virtual void centroid( const Teuchos::ArrayView<double>& centroid ) const;

    /*!
     * \brief Return the axis-aligned bounding box around the entity.
     * \param bounding_box A Cartesian box that bounds the entity.
     */
    virtual void boundingBox( Box& bounding_box ) const;
    //@}

    //@{
    //! Parameteric mapping functions.
    /*!
     * \brief Perform a safeguard check for mapping a point to the reference
     * space of an entity using the given tolerance.
     * \param parameters Parameters to be used for the safeguard check.
     * \param point A view into an array of size physicalDimension() containing
     * the coordinates of the point to map.
     * \param status A status object indicating the results of the safeguard
     * check.
     */
    virtual void safeguardMapToReferenceFrame(
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& point,
	MappingStatus& status ) const;

    /*!
     * \brief Map a point to the reference space of an entity. Return the
     * parameterized point.
     * \param parameters Parameters to be used for the mapping procedure.
     * \param  A view into an array of size physicalDimension() containing
     * the coordinates of the point to map.
     * \param reference_point A view into an array of size physicalDimension()
     * to write the reference coordinates of the mapped point.
     * \param status A status object indicating the results of the mapping
     * procedure.
     */
    virtual void mapToReferenceFrame( 
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point,
	MappingStatus& status ) const;

    /*!  
     * \brief Determine if a reference point is in the parameterized space of
     * an entity.
     * \param parameters Parameters to be used for the point inclusion check.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \return True if the point is in the reference space, false if not.
    */
    virtual bool checkPointInclusion( 
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& reference_point ) const;

    /*!
     * \brief Map a reference point to the physical space of an entity.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \param  A view into an array of size physicalDimension() to write
     * the coordinates of physical point.
     */
    virtual void mapToPhysicalFrame( 
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point ) const;
    //@}

    //@{
    //! SerializableAbstractObjectPolicy interface.
    /*!
     * \brief Return a string indicating the derived object type.
     * \return A string indicating the type of derived object implementing the
     * interface.
     */
    std::string objectType() const;

    /*
     * \brief Get the size of the serialized entity in bytes.
     * \return The size of the entity when serialized in bytes.
     */
    virtual std::size_t byteSize() const;

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

    //@{
    //! Static polymorphic construction interface.
    /*!
     * \brief Set an abstract factor for a GeometricEntity subclass.
     * \param builder A factory for a GeometricEntity subclass.
     */
    static void setDerivedClassFactory(
	const Teuchos::RCP<
	const Teuchos::AbstractFactory<GeometricEntity> >& factory );

    /*!
     * \brief Get the abstract builder for GeometricEntity subclasses.
     * \return The builder for GeometricEntity subclasses.
     */
    static Teuchos::RCP<AbstractBuilder<GeometricEntity> > getBuilder();
    //@}

  private:
    
    // Abstract builder for geometric entity subclasses.
    static Teuchos::RCP<AbstractBuilder<GeometricEntity> > b_builder;
};

//---------------------------------------------------------------------------//
// SerializableAbstractObjectPolicy implementation.
//---------------------------------------------------------------------------//

template<>
class SerializableAbstractObjectPolicy<GeometricEntity>
{
  public:

    typedef GeometricEntity object_type;

    static std::string objectType( const Teuchos::RCP<object_type>& object )
    {
	return object->objectType();
    }

    static std::size_t byteSize( const Teuchos::RCP<object_type>& object )
    {
	return object->byteSize();
    }

    static void serialize( const Teuchos::RCP<object_type>& object,
			   const Teuchos::ArrayView<char>& buffer )
    {
	object->serialize( buffer );
    }

    static void deserialize( const Teuchos::RCP<object_type>& object,
			     const Teuchos::ArrayView<const char>& buffer )
    {
	object->deserialize( buffer );
    }

    static Teuchos::RCP<DataTransferKit::AbstractBuilder<object_type> > getBuilder()
    {
	return object_type::getBuilder();
    }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Teuchos::SerializationTraits implementation.
//---------------------------------------------------------------------------//

namespace Teuchos
{
template<typename Ordinal>
class SerializationTraits<Ordinal,Teuchos::RCP<DataTransferKit::GeometricEntity> > 
{
  public:

    typedef DataTransferKit::GeometricEntity Base;
    typedef Teuchos::RCP<Base> T;
    typedef DataTransferKit::AbstractSerializer<Ordinal,Base>  
    AbstractSerializer;

    static const bool supportsDirectSerialization = 
	AbstractSerializer::supportsDirectSerialization;

    static Ordinal fromCountToIndirectBytes( const Ordinal count, 
					     const T buffer[] ) 
    { 
	return AbstractSerializer::fromCountToIndirectBytes( count, buffer );
    }

    static void serialize( const Ordinal count, 
			   const T buffer[], 
			   const Ordinal bytes, 
			   char charBuffer[] )
    { 
	AbstractSerializer::serialize( count, buffer, bytes, charBuffer );
    }

    static Ordinal fromIndirectBytesToCount( const Ordinal bytes, 
					     const char charBuffer[] ) 
    { 
	return AbstractSerializer::fromIndirectBytesToCount( bytes, charBuffer );
    }

    static void deserialize( const Ordinal bytes, 
			     const char charBuffer[], 
			     const Ordinal count, 
			     T buffer[] )
    { 
	AbstractSerializer::deserialize( bytes, charBuffer, count, buffer );
    }
};
} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end DTK_GEOMETRICENTITY_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometricEntity.hpp
//---------------------------------------------------------------------------//
