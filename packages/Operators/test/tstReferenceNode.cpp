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
 * \file tstReferenceNode.cpp
 * \author Stuart R. Slattery
 * \brief ReferenceNode unit tests.
 */
//---------------------------------------------------------------------------//

#include "reference_implementation/DTK_ReferenceNode.hpp"

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ReferenceNode, reference_node )
{
    DataTransferKit::EntityId id = 3242;
    int owner_rank = 16;
    double x = -33.2;
    double y = 147.6;
    double z = 55.9;

    DataTransferKit::Entity entity =
        DataTransferKit::UnitTest::ReferenceNode( id, owner_rank, x, y, z );

    TEST_EQUALITY( entity.id(), id );
    TEST_EQUALITY( entity.ownerRank(), owner_rank );
    TEST_EQUALITY( entity.topologicalDimension(), 0 );
    TEST_EQUALITY( entity.physicalDimension(), 3 );

    Teuchos::Tuple<double,6> box;
    entity.boundingBox( box );
    TEST_EQUALITY( box[0], x );
    TEST_EQUALITY( box[1], y );
    TEST_EQUALITY( box[2], z );
    TEST_EQUALITY( box[3], x );
    TEST_EQUALITY( box[4], y );
    TEST_EQUALITY( box[5], z );

    TEST_ASSERT( !entity.inBlock(0) );
    TEST_ASSERT( !entity.onBoundary(0) );

    std::cout << entity.description() << std::endl;

    Teuchos::RCP<Teuchos::FancyOStream>
        fancy_out = Teuchos::VerboseObjectBase::getDefaultOStream();
    entity.describe( *fancy_out );
}

//---------------------------------------------------------------------------//
// end tstReferenceNode.cpp
//---------------------------------------------------------------------------//

