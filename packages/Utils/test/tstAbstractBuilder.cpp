//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstAbstractBuilder.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  AbstractBuilder class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DTK_AbstractBuilder.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <Teuchos_AbstractFactoryStd.hpp>

//---------------------------------------------------------------------------//
// HELPER CLASSES
//---------------------------------------------------------------------------//

class BaseClass
{
  public:
    BaseClass() { /* ... */ }
    virtual ~BaseClass() { /* ... */ }
    virtual int myNumber() = 0;
};

class MyNumberIsOne : public BaseClass
{
  public:
    MyNumberIsOne() { /* ... */ }
    ~MyNumberIsOne() { /* ... */ }
    int myNumber() { return 1; }
};

class MyNumberIsTwo : public BaseClass
{
  public:
    MyNumberIsTwo() { /* ... */ }
    ~MyNumberIsTwo() { /* ... */ }
    int myNumber() { return 2; }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( AbstractBuilder, builder_test )
{
    DataTransferKit::AbstractBuilder<BaseClass> builder;

    std::string name_1( "one" );
    builder.setDerivedClassFactory(
	Teuchos::abstractFactoryStd<BaseClass,MyNumberIsOne>(), name_1 );

    std::string name_2( "two" );
    builder.setDerivedClassFactory(
	Teuchos::abstractFactoryStd<BaseClass,MyNumberIsTwo>(), name_2 );

    Teuchos::RCP<BaseClass> base_1 = builder.create( "one" );
    TEST_EQUALITY( base_1->myNumber(), 1 );

    Teuchos::RCP<BaseClass> base_2 = builder.create( "two" );
    TEST_EQUALITY( base_2->myNumber(), 2 );
}

//---------------------------------------------------------------------------//
//                        end of tstAbstractBuilder.cpp
//---------------------------------------------------------------------------//
