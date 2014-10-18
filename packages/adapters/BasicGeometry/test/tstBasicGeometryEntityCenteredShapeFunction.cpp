//---------------------------------------------------------------------------//
/*!
 * \file tstBasicGeometryEntityCenteredShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief BasicGeometryEntityCenteredShapeFunction unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_EntityShapeFunction.hpp>
#include <DTK_BasicGeometryEntityCenteredShapeFunction.hpp>
#include <DTK_Point.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( BasicGeometryEntityCenteredShapeFunction, shape_func_test )
{
    using namespace DataTransferKit;

    // Make point.
    Teuchos::Array<double> p(3);
    p[0] = 3.2;
    p[1] = 302.3;
    p[2] = 9.32;
    EntityId id = 12;
    Entity point = Point( id, 0, p);

    // Make a shape function.
    Teuchos::RCP<EntityShapeFunction> shape_function =
	Teuchos::rcp( new BasicGeometryEntityCenteredShapeFunction() );

    // Test the shape function.
    Teuchos::Array<std::size_t> dof_ids;
    shape_function->entityDOFIds( point, dof_ids );
    TEST_EQUALITY( 1, dof_ids.size() );
    TEST_EQUALITY( Teuchos::as<std::size_t>(id), dof_ids[0] );

    Teuchos::Array<double> values;
    shape_function->evaluateValue( point, p(), values );
    TEST_EQUALITY( 1, values.size() );
    TEST_EQUALITY( 1.0, values[0] );

    Teuchos::Array<Teuchos::Array<double> > gradients;
    shape_function->evaluateGradient( point, p(), gradients );
    TEST_EQUALITY( 1, gradients.size() );
    TEST_EQUALITY( 1, gradients[0].size() );
    TEST_EQUALITY( 0.0, gradients[0][0] );
}

//---------------------------------------------------------------------------//
// end tstBasicGeometryEntityCenteredShapeFunction.cpp
//---------------------------------------------------------------------------//

