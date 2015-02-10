//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file tstEntityCenteredShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief EntityCenteredShapeFunction unit tests.
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
#include <DTK_EntityCenteredShapeFunction.hpp>
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
TEUCHOS_UNIT_TEST( EntityCenteredShapeFunction, shape_func_test )
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
	Teuchos::rcp( new EntityCenteredShapeFunction() );

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
// end tstEntityCenteredShapeFunction.cpp
//---------------------------------------------------------------------------//

