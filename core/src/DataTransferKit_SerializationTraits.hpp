//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DataTransferKit_SerializationTraits.hpp
 * \author Stuart Slattery
 * \date   Fri Dec 09 09:25:58 2011
 * \brief  Serialization traits for datatransferkit objects.
 */
//---------------------------------------------------------------------------//

#ifndef DATATRANSFERKIT_SERIALIZATIONTRAITS_HPP
#define DATATRANSFERKIT_SERIALIZATIONTRAITS_HPP

#include "DataTransferKit_Point.hpp"
#include "DataTransferKit_BoundingBox.hpp"

#include <Teuchos_SerializationTraits.hpp>

namespace Teuchos 
{
//===========================================================================//
/*!
 * \brief Serialization traits for datatransferkit objects.
 */
//===========================================================================//

// Point.
template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<long int,float> >
{};

// Bounding box.
template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::BoundingBox<int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::BoundingBox<int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::BoundingBox<int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::BoundingBox<int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::BoundingBox<long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::BoundingBox<long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::BoundingBox<long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::BoundingBox<long int,float> >
{};

} // end namepsace Teuchos

#endif // DATATRANSFERKIT_SERIALIZATIONTRAITS_HPP

//---------------------------------------------------------------------------//
//              end of DataTransferKit_SerializationTraits.hpp
//---------------------------------------------------------------------------//
