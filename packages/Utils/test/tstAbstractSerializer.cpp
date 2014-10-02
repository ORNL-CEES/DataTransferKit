//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstAbstractSerializer.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Abstract serializer class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DTK_AbstractSerializer.hpp>
#include <DTK_AbstractSerializableObjectPolicy.hpp>
#include <DTK_AbstractBuilder.hpp>\

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_SerializationTraits.hpp"

//---------------------------------------------------------------------------//
// HELPER CLASSES
//---------------------------------------------------------------------------//

// Base class.
class BaseClass
{
  public:

    BaseClass() { /* ... */ }
    virtual ~BaseClass() { /* ... */ }

    virtual int myNumber() { return 0; }
    virtual Teuchos::Array<double> myData() { return Teuchos::Array<double>(0); }
        
    virtual std::string objectType() const { return std::string("0"); }
    virtual std::size_t byteSize() const { return 0; }
    virtual void serialize( const Teuchos::ArrayView<char>& buffer ) const {};
    virtual void deserialize( const Teuchos::ArrayView<const char>& buffer ) {};

    static void setBuilder( 
	const Teuchos::RCP<DataTransferKit::AbstractBuilder<BaseClass> >& builder );
    static Teuchos::RCP<DataTransferKit::AbstractBuilder<BaseClass> > getBuilder();

  private:
    
    static Teuchos::RCP<DataTransferKit::AbstractBuilder<BaseClass> > d_builder;
};

Teuchos::RCP<DataTransferKit::AbstractBuilder<BaseClass> > 
BaseClass::d_builder = Teuchos::null;

void BaseClass::setBuilder( 
    const Teuchos::RCP<DataTransferKit::AbstractBuilder<BaseClass> >& builder )
{ d_builder = builder; }

Teuchos::RCP<DataTransferKit::AbstractBuilder<BaseClass> > BaseClass::getBuilder()
{ return d_builder; }

//---------------------------------------------------------------------------//
// AbstractSerializableObjectPolicy implementation for the base class.
namespace DataTransferKit
{
template<>
class AbstractSerializableObjectPolicy<BaseClass>
{
  public:

    typedef BaseClass object_type;

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
// SerializationTraits implementation for the base class pointers.
namespace Teuchos
{
template<typename Ordinal>
class SerializationTraits<Ordinal,Teuchos::RCP<BaseClass> > 
{
  public:

    typedef Teuchos::RCP<BaseClass> T;
    typedef DataTransferKit::AbstractSerializer<Ordinal,BaseClass>  
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
// Derived class 1
class MyNumberIsOne : public BaseClass
{
  public:

    MyNumberIsOne() : d_data( 1, 1.0 ) { /* ... */ }
    ~MyNumberIsOne() { /* ... */ }

    int myNumber() { return 1; }
    Teuchos::Array<double> myData() { return d_data; }

    std::string objectType() const { return std::string("one"); }
    std::size_t byteSize() const { return sizeof(double); }
    void serialize( const Teuchos::ArrayView<char>& buffer ) const
    { std::memcpy( const_cast<double*>(d_data.getRawPtr()), 
		   buffer.getRawPtr(), sizeof(double) ); }
    void deserialize( const Teuchos::ArrayView<const char>& buffer )
    { std::memcpy( const_cast<char*>(buffer.getRawPtr()), 
		   d_data.getRawPtr(), sizeof(double) ); }

  private: 
    
    Teuchos::Array<double> d_data;
};

//---------------------------------------------------------------------------//
// Derived class 2.
class MyNumberIsTwo : public BaseClass
{
  public:

    MyNumberIsTwo() : d_data( 2, 2.0 ) { /* ... */ }
    ~MyNumberIsTwo() { /* ... */ }

    int myNumber() { return 2; }
    Teuchos::Array<double> myData() { return d_data; }

    std::string objectType() const { return std::string("two"); }
    std::size_t byteSize() const { return 2*sizeof(double); }
    void serialize( const Teuchos::ArrayView<char>& buffer ) const
    { std::memcpy( const_cast<double*>(d_data.getRawPtr()), 
		   buffer.getRawPtr(), sizeof(double) ); }
    void deserialize( const Teuchos::ArrayView<const char>& buffer )
    { std::memcpy( const_cast<char*>(buffer.getRawPtr()), 
		   d_data.getRawPtr(), sizeof(double) ); }

  private: 
    
    Teuchos::Array<double> d_data;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( AbstractSerializableObjectPolicy, serializable_object_test )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_default = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm_default->getRank();

    // Create a builder and register derived classes.
    Teuchos::RCP<AbstractBuilder<BaseClass> > builder =
	Teuchos::rcp( new AbstractBuilder<BaseClass>() );

    std::string name_1( "one" );
    builder->setDerivedClassFactory(
	Teuchos::abstractFactoryStd<BaseClass,MyNumberIsOne>(), name_1 );

    std::string name_2( "two" );
    builder->setDerivedClassFactory(
	Teuchos::abstractFactoryStd<BaseClass,MyNumberIsTwo>(), name_2 );

    // Set the builder with the base class.
    BaseClass::setBuilder( builder );

    // Construct an array of base class objects.
    Teuchos::Array<Teuchos::RCP<BaseClass> > objects( 2 );
    if ( comm_rank == 0 )
    {
	objects[0] = builder->create("one");
	objects[1] = builder->create("two");
    }

    // Broadcast the objects.
    Teuchos::broadcast( *comm_default, 0, objects() );

    // Check the objects.
    TEST_EQUALITY( 1, objects[0]->myNumber() );
    TEST_EQUALITY( 1, objects[0]->myData().size() );
    TEST_EQUALITY( 1.0, objects[0]->myData()[0] );

    TEST_EQUALITY( 2, objects[1]->myNumber() );
    TEST_EQUALITY( 2, objects[1]->myData().size() );
    TEST_EQUALITY( 2.0, objects[1]->myData()[0] );
    TEST_EQUALITY( 2.0, objects[1]->myData()[1] );
}

//---------------------------------------------------------------------------//
// end of tstAbstractSerializer.cpp
//---------------------------------------------------------------------------//
