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
#include <DTK_BasicEntitySet.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>
#include <DTK_EntityCenteredShapeFunction.hpp>

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

    Teuchos::RCP<EntitySet> entity_set = Teuchos::rcp(
	new BasicEntitySet(Teuchos::DefaultComm<int>::getComm(), 1) );
    Teuchos::RCP<EntityLocalMap> local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );
    Teuchos::RCP<EntityShapeFunction> shape_function =
	Teuchos::rcp( new EntityCenteredShapeFunction() );
    FunctionSpace function_space( entity_set, local_map, shape_function );

    TEST_EQUALITY( function_space.entitySet().getRawPtr(), 
		   entity_set.getRawPtr() );
    TEST_EQUALITY( function_space.localMap().getRawPtr(),
		   local_map.getRawPtr() );
    TEST_EQUALITY( function_space.shapeFunction().getRawPtr(),
		   shape_function.getRawPtr() );
}

//---------------------------------------------------------------------------//
// end tstFunctionSpace.cpp
//---------------------------------------------------------------------------//
