//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Coupler_SerializationTraits.hpp
 * \author Stuart Slattery
 * \date   Fri Dec 09 09:25:58 2011
 * \brief  Serialization traits for coupler objects.
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_SERIALIZATIONTRAITS_HPP
#define COUPLER_SERIALIZATIONTRAITS_HPP

#include "Coupler_Point.hpp"

#include <Teuchos_SerializationTraits.hpp>

namespace Teuchos 
{
//===========================================================================//
/*!
 * \brief Serialization traits for coupler objects.
 */
//===========================================================================//

// Point.
template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<1,int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<1,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<1,int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<1,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<1,long int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<1,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<1,long int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<1,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<2,int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<2,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<2,int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<2,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<2,long int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<2,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<2,long int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<2,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<3,int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<3,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<3,int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<3,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<3,long int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<3,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<3,long int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<3,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<4,int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<4,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<4,int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<4,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<4,long int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<4,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<4,long int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<4,long int,float> >
{};

} // end namepsace Teuchos

#endif // COUPLER_SERIALIZATIONTRAITS_HPP

//---------------------------------------------------------------------------//
//              end of Coupler_SerializationTraits.hpp
//---------------------------------------------------------------------------//
