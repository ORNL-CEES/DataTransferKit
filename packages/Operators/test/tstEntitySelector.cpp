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
#include <DTK_EntityPredicates.hpp>
#include <DTK_PredicateComposition.hpp>
#include <DTK_Point.hpp>

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

    Teuchos::Array<double> coords1(1,0.0);
    Teuchos::Array<int> blocks1(0);
    Teuchos::Array<int> boundaries1(0);
    Entity p1 = Point( 0, 0, coords1, false, blocks1, boundaries1 );

    Teuchos::Array<double> coords2(2,0.0);
    Teuchos::Array<int> blocks2(0);
    Teuchos::Array<int> boundaries2(0);
    Entity p2 = Point( 0, 0, coords2, true, blocks2, boundaries2 );

    EntitySelector select1( ENTITY_TYPE_NODE );
    TEST_EQUALITY( ENTITY_TYPE_NODE, select1.entityType() );
    TEST_ASSERT( select1.selectFunction()(p1) );
    TEST_ASSERT( select1.selectFunction()(p2) );

    SurfacePredicate surf_pred;
    EntitySelector select2( ENTITY_TYPE_NODE, surf_pred );
    TEST_EQUALITY( ENTITY_TYPE_NODE, select2.entityType() );
    TEST_ASSERT( !select2.selectFunction()(p1) );
    TEST_ASSERT( select2.selectFunction()(p2) );
}

//---------------------------------------------------------------------------//
// end tstEntitySelector.cpp
//---------------------------------------------------------------------------//
