//---------------------------------------------------------------------------//
/*! 
 * \file tstFunctionSpace.cpp
 * \author Stuart R. Slattery
 * \brief FunctionSpace unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_FunctionSpace.hpp>
#include <DTK_EntitySet.hpp>
#include <DTK_EntityLocalMap.hpp>
#include <DTK_EntityShapeFunction.hpp>
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
TEUCHOS_UNIT_TEST( FunctionSpace, space_test )
{
    using namespace DataTransferKit;

    Teuchos::RCP<EntitySet> entity_set = Teuchos::rcp( new EntitySet() );
    Teuchos::RCP<EntityLocalMap> local_map = Teuchos::rcp( new EntityLocalMap() );
    Teuchos::RCP<EntityShapeFunction> shape_function =
	Teuchos::rcp( new EntityShapeFunction() );
    Teuchos::RCP<EntitySelector> selector =
	Teuchos::rcp( new EntitySelector(ENTITY_TYPE_NODE) );
    FunctionSpace function_space( entity_set, selector, local_map, shape_function );

    TEST_EQUALITY( function_space.entitySet().getRawPtr(), 
		   entity_set.getRawPtr() );
    TEST_EQUALITY( function_space.entitySelector().getRawPtr(), 
		   selector.getRawPtr() );
    TEST_EQUALITY( function_space.localMap().getRawPtr(),
		   local_map.getRawPtr() );
    TEST_EQUALITY( function_space.shapeFunction().getRawPtr(),
		   shape_function.getRawPtr() );
}

//---------------------------------------------------------------------------//
// end tstFunctionSpace.cpp
//---------------------------------------------------------------------------//
