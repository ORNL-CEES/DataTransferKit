//---------------------------------------------------------------------------//
/*! 
 * \file tstBasicEntityPredicates.cpp
 * \author Stuart R. Slattery
 * \brief BasicEntityPredicates unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_BasicEntityPredicates.hpp>
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
TEUCHOS_UNIT_TEST( BlockPredicate, block_predicate_test )
{
    using namespace DataTransferKit;

    Teuchos::Array<double> coords1(1,0.0);
    Teuchos::Array<int> blocks1(1,2);
    Teuchos::Array<int> boundaries1(0);
    Entity p1 = Point( 0, 0, coords1, blocks1, boundaries1 );

    Teuchos::Array<double> coords2(2,0.0);
    Teuchos::Array<int> blocks2(1,1);
    Teuchos::Array<int> boundaries2(0);
    Entity p2 = Point( 0, 0, coords2, blocks2, boundaries2 );

    Teuchos::Array<int> pred_1(1,1);
    BlockPredicate block_pred_1( pred_1 );
    TEST_ASSERT( !block_pred_1(p1) );
    TEST_ASSERT( block_pred_1(p2) );

    Teuchos::Array<int> pred_2(1,2);
    BlockPredicate block_pred_2( pred_2 );
    TEST_ASSERT( block_pred_2(p1) );
    TEST_ASSERT( !block_pred_2(p2) );

    Teuchos::Array<int> pred_3(2);
    pred_3[0] = 1;
    pred_3[1] = 2;
    BlockPredicate block_pred_3( pred_3 );
    TEST_ASSERT( block_pred_3(p1) );
    TEST_ASSERT( block_pred_3(p2) );

    std::function<bool(Entity)> block_pred_4 =
	PredicateComposition::And( block_pred_1.getFunction(), 
				   block_pred_2.getFunction() );
    TEST_ASSERT( !block_pred_4(p1) );
    TEST_ASSERT( !block_pred_4(p2) );

    std::function<bool(Entity)> block_pred_5 =
	PredicateComposition::Or( block_pred_1.getFunction(), 
				  block_pred_2.getFunction() );
    TEST_ASSERT( block_pred_5(p1) );
    TEST_ASSERT( block_pred_5(p2) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( BoundPredicate, bound_predicate_test )
{
    using namespace DataTransferKit;

    Teuchos::Array<double> coords1(1,0.0);
    Teuchos::Array<int> bounds1(0);
    Teuchos::Array<int> boundaries1(1,2);
    Entity p1 = Point( 0, 0, coords1, bounds1, boundaries1 );

    Teuchos::Array<double> coords2(2,0.0);
    Teuchos::Array<int> bounds2(0);
    Teuchos::Array<int> boundaries2(1,1);
    Entity p2 = Point( 0, 0, coords2, bounds2, boundaries2 );

    Teuchos::Array<int> pred_1(1,1);
    BoundaryPredicate bound_pred_1( pred_1 );
    TEST_ASSERT( !bound_pred_1(p1) );
    TEST_ASSERT( bound_pred_1(p2) );

    Teuchos::Array<int> pred_2(1,2);
    BoundaryPredicate bound_pred_2( pred_2 );
    TEST_ASSERT( bound_pred_2(p1) );
    TEST_ASSERT( !bound_pred_2(p2) );

    Teuchos::Array<int> pred_3(2);
    pred_3[0] = 1;
    pred_3[1] = 2;
    BoundaryPredicate bound_pred_3( pred_3 );
    TEST_ASSERT( bound_pred_3(p1) );
    TEST_ASSERT( bound_pred_3(p2) );

    std::function<bool(Entity)> bound_pred_4 =
	PredicateComposition::And( bound_pred_1.getFunction(),
				   bound_pred_2.getFunction() );
    TEST_ASSERT( !bound_pred_4(p1) );
    TEST_ASSERT( !bound_pred_4(p2) );

    std::function<bool(Entity)> bound_pred_5 =
	PredicateComposition::Or( bound_pred_1.getFunction(), 
				  bound_pred_2.getFunction() );
    TEST_ASSERT( bound_pred_5(p1) );
    TEST_ASSERT( bound_pred_5(p2) );
}

//---------------------------------------------------------------------------//
// end tstBasicEntityPredicates.cpp
//---------------------------------------------------------------------------//
