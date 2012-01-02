//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Mesh_SerializationTraits.hpp
 * \author Stuart Slattery
 * \date   Fri Dec 09 09:25:58 2011
 * \brief  Serialization traits for mesh objects.
 */
//---------------------------------------------------------------------------//

#ifndef mesh_Mesh_SerializationTraits_hpp
#define mesh_Mesh_SerializationTraits_hpp

#include <Mesh_Point.hpp>
#include <Mesh_BoundingBox.hpp>

#include "Teuchos_SerializationTraits.hpp"

//===========================================================================//
/*!
 * \class Mesh_SerializationTraits
 * \brief Serialization traits for mesh objects.
 */
//===========================================================================//

namespace Teuchos {

// Point.
template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<long int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::Point<long int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::Point<long int,float> >
{};

// Bounding box.
template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::BoundingBox<int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::BoundingBox<int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::BoundingBox<int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::BoundingBox<int,float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::BoundingBox<long int,double> >
    : public DirectSerializationTraits<Ordinal,Coupler::BoundingBox<long int,double> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,Coupler::BoundingBox<long int,float> >
    : public DirectSerializationTraits<Ordinal,Coupler::BoundingBox<long int,float> >
{};

} // end namepsace Teuchos

#endif // mesh_Mesh_SerializationTraits_hpp

//---------------------------------------------------------------------------//
//              end of mesh/Mesh_SerializationTraits.hpp
//---------------------------------------------------------------------------//
