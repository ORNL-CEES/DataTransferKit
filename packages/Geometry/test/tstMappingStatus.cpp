//---------------------------------------------------------------------------//
/*!
 * \file tstMappingStatus.cpp
 * \author Stuart R. Slattery
 * \brief Bounding Box unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MappingStatus.hpp>

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
TEUCHOS_UNIT_TEST( Box, status_test )
{
    DataTransferKit::MappingStatus status_1;
    TEST_EQUALITY( status_1.success(), false );
    TEST_EQUALITY( status_1.numberOfIterations(), 0 );

    status_1.mappingSucceeded();
    TEST_EQUALITY( status_1.success(), true );

    status_1.mappingFailed();
    TEST_EQUALITY( status_1.success(), false );

    unsigned iterations = 100;
    for ( unsigned i = 0; i < iterations; ++i )
    {
	status_1.incrementIterations();
    }
    TEST_EQUALITY( status_1.numberOfIterations(), iterations );

    for ( unsigned i = 0; i < iterations; ++i )
    {
	status_1.incrementIterations();
    }
    TEST_EQUALITY( status_1.numberOfIterations(), 2*iterations );

    status_1.reset();
    TEST_EQUALITY( status_1.success(), false );
    TEST_EQUALITY( status_1.numberOfIterations(), 0 );

    DataTransferKit::MappingStatus status_2( true, iterations );
    TEST_EQUALITY( status_2.success(), true );
    TEST_EQUALITY( status_2.numberOfIterations(), iterations );

    DataTransferKit::MappingStatus status_3 = status_2;
    TEST_EQUALITY( status_3.success(), true );
    TEST_EQUALITY( status_3.numberOfIterations(), iterations );
}

//---------------------------------------------------------------------------//
// end tstMappingStatus.cpp
//---------------------------------------------------------------------------//

