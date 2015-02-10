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
