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
 * \file tstBasicEntityPredicates.cpp
 * \author Stuart R. Slattery
 * \brief BasicEntityPredicates unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_BasicEntityPredicates.hpp>
#include <DTK_Entity.hpp>
#include <DTK_EntityImpl.hpp>
#include <DTK_PredicateComposition.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// Basic entity implementation
//---------------------------------------------------------------------------//
class MyEntityImpl : public DataTransferKit::EntityImpl
{
  public:
    MyEntityImpl() { /* ... */}
    MyEntityImpl( const Teuchos::Array<int> &block_ids,
                  const Teuchos::Array<int> &boundary_ids )
        : d_block_ids( block_ids )
        , d_boundary_ids( boundary_ids )
    {
        std::sort( d_block_ids.begin(), d_block_ids.end() );
        std::sort( d_boundary_ids.begin(), d_boundary_ids.end() );
    }
    DataTransferKit::EntityId id() const { return 0; }
    int ownerRank() const { return 0; }
    int topologicalDimension() const { return 0; }
    int physicalDimension() const { return 0; }
    void boundingBox( Teuchos::Tuple<double, 6> &bounds ) const {}
    bool inBlock( const int block_id ) const
    {
        return std::binary_search( d_block_ids.begin(), d_block_ids.end(),
                                   block_id );
    }
    bool onBoundary( const int boundary_id ) const
    {
        return std::binary_search( d_boundary_ids.begin(), d_boundary_ids.end(),
                                   boundary_id );
    }
    Teuchos::RCP<DataTransferKit::EntityExtraData> extraData() const
    {
        return Teuchos::null;
    }

  private:
    Teuchos::Array<int> d_block_ids;
    Teuchos::Array<int> d_boundary_ids;
};

class MyEntity : public DataTransferKit::Entity
{
  public:
    MyEntity() { /* ... */}
    MyEntity(
        const Teuchos::Array<int> &block_ids = Teuchos::Array<int>( 0 ),
        const Teuchos::Array<int> &boundary_ids = Teuchos::Array<int>( 0 ) )
    {
        this->b_entity_impl =
            Teuchos::rcp( new MyEntityImpl( block_ids, boundary_ids ) );
    }
};

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( BlockPredicate, block_predicate_test )
{
    using namespace DataTransferKit;

    Teuchos::Array<int> blocks1( 1, 2 );
    Teuchos::Array<int> boundaries1( 0 );
    Entity p1 = MyEntity( blocks1, boundaries1 );

    Teuchos::Array<int> blocks2( 1, 1 );
    Teuchos::Array<int> boundaries2( 0 );
    Entity p2 = MyEntity( blocks2, boundaries2 );

    Teuchos::Array<int> pred_1( 1, 1 );
    BlockPredicate block_pred_1( pred_1 );
    TEST_ASSERT( !block_pred_1( p1 ) );
    TEST_ASSERT( block_pred_1( p2 ) );

    Teuchos::Array<int> pred_2( 1, 2 );
    BlockPredicate block_pred_2( pred_2 );
    TEST_ASSERT( block_pred_2( p1 ) );
    TEST_ASSERT( !block_pred_2( p2 ) );

    Teuchos::Array<int> pred_3( 2 );
    pred_3[0] = 1;
    pred_3[1] = 2;
    BlockPredicate block_pred_3( pred_3 );
    TEST_ASSERT( block_pred_3( p1 ) );
    TEST_ASSERT( block_pred_3( p2 ) );

    std::function<bool( Entity )> block_pred_4 = PredicateComposition::And(
        block_pred_1.getFunction(), block_pred_2.getFunction() );
    TEST_ASSERT( !block_pred_4( p1 ) );
    TEST_ASSERT( !block_pred_4( p2 ) );

    std::function<bool( Entity )> block_pred_5 = PredicateComposition::Or(
        block_pred_1.getFunction(), block_pred_2.getFunction() );
    TEST_ASSERT( block_pred_5( p1 ) );
    TEST_ASSERT( block_pred_5( p2 ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( BoundPredicate, bound_predicate_test )
{
    using namespace DataTransferKit;

    Teuchos::Array<int> bounds1( 0 );
    Teuchos::Array<int> boundaries1( 1, 2 );
    Entity p1 = MyEntity( bounds1, boundaries1 );

    Teuchos::Array<int> bounds2( 0 );
    Teuchos::Array<int> boundaries2( 1, 1 );
    Entity p2 = MyEntity( bounds2, boundaries2 );

    Teuchos::Array<int> pred_1( 1, 1 );
    BoundaryPredicate bound_pred_1( pred_1 );
    TEST_ASSERT( !bound_pred_1( p1 ) );
    TEST_ASSERT( bound_pred_1( p2 ) );

    Teuchos::Array<int> pred_2( 1, 2 );
    BoundaryPredicate bound_pred_2( pred_2 );
    TEST_ASSERT( bound_pred_2( p1 ) );
    TEST_ASSERT( !bound_pred_2( p2 ) );

    Teuchos::Array<int> pred_3( 2 );
    pred_3[0] = 1;
    pred_3[1] = 2;
    BoundaryPredicate bound_pred_3( pred_3 );
    TEST_ASSERT( bound_pred_3( p1 ) );
    TEST_ASSERT( bound_pred_3( p2 ) );

    std::function<bool( Entity )> bound_pred_4 = PredicateComposition::And(
        bound_pred_1.getFunction(), bound_pred_2.getFunction() );
    TEST_ASSERT( !bound_pred_4( p1 ) );
    TEST_ASSERT( !bound_pred_4( p2 ) );

    std::function<bool( Entity )> bound_pred_5 = PredicateComposition::Or(
        bound_pred_1.getFunction(), bound_pred_2.getFunction() );
    TEST_ASSERT( bound_pred_5( p1 ) );
    TEST_ASSERT( bound_pred_5( p2 ) );
}

//---------------------------------------------------------------------------//
// end tstBasicEntityPredicates.cpp
//---------------------------------------------------------------------------//
