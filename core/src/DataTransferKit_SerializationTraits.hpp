//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DataTransferKit_SerializationTraits.hpp
 * \author Stuart Slattery
 * \date   Fri Dec 09 09:25:58 2011
 * \brief  Serialization traits for core objects.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SERIALIZATIONTRAITS_HPP
#define DTK_SERIALIZATIONTRAITS_HPP

#include "DataTransferKit_Point.hpp"

#include <Teuchos_SerializationTraits.hpp>

namespace Teuchos 
{
//===========================================================================//
/*!
 * \brief Serialization traits for core objects.
 */
//===========================================================================//

// Point.
template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<1,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<1,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<1,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<1,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<1,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<1,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<1,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<1,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<2,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<2,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<2,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<2,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<2,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<2,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<2,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<2,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<3,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<3,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<3,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<3,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<3,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<3,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<3,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<3,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<4,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<4,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<4,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<4,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<4,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<4,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Point<4,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Point<4,long int,float> >
{};

} // end namepsace Teuchos

#endif // DTK_SERIALIZATIONTRAITS_HPP

//---------------------------------------------------------------------------//
//              end of DataTransferKit_SerializationTraits.hpp
//---------------------------------------------------------------------------//
