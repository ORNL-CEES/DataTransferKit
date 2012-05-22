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

#include "DTK_Node.hpp"
#include "DTK_Element.hpp"

#include <Teuchos_SerializationTraits.hpp>

namespace Teuchos 
{
//===========================================================================//
/*!
 * \brief Serialization traits for core objects.
 */
//===========================================================================//

// Node.
template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<1,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<1,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<1,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<1,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<1,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<1,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<1,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<1,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<2,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<2,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<2,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<2,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<2,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<2,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<2,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<2,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<3,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<3,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<3,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<3,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<3,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<3,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<3,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<3,long int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<4,int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<4,int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<4,int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<4,int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<4,long int,double> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<4,long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Node<4,long int,float> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Node<4,long int,float> >
{};

// Element.
template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<1,int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<1,int> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<1,long int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<1,long int> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<2,int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<2,int> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<2,long int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<2,long int> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<3,int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<3,int> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<3,long int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<3,long int> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<4,int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<4,int> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Element<4,long int> >
    : public DirectSerializationTraits<Ordinal,DataTransferKit::Element<4,long int> >
{};

} // end namepsace Teuchos

#endif // DTK_SERIALIZATIONTRAITS_HPP

//---------------------------------------------------------------------------//
//              end of DataTransferKit_SerializationTraits.hpp
//---------------------------------------------------------------------------//
