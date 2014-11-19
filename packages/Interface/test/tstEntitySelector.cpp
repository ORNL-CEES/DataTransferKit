//---------------------------------------------------------------------------//
/*! 
 * \file tstEntitySelector.cpp
 * \author Stuart R. Slattery
 * \brief EntitySelector unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_EntitySelector.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( EntitySelector, selector_test )
{
    using namespace DataTransferKit;

    Entity p1;
    Entity p2;

    EntitySelector select1( ENTITY_TYPE_NODE );
    TEST_EQUALITY( ENTITY_TYPE_NODE, select1.entityType() );
    TEST_ASSERT( select1.selectFunction()(p1) );
    TEST_ASSERT( select1.selectFunction()(p2) );
}

//---------------------------------------------------------------------------//
// end tstEntitySelector.cpp
//---------------------------------------------------------------------------//
