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
 * \file tstEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief Abstract iterator unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_Entity.hpp>
#include <DTK_EntityIterator.hpp>
#include <DTK_PredicateComposition.hpp>
#include <DTK_Types.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// EntityImpl Implementation.
//---------------------------------------------------------------------------//
class TestEntityImpl : public DataTransferKit::EntityImpl
{
  public:
    TestEntityImpl( int id ) { d_id = id; }
    DataTransferKit::EntityId id() const { return d_id; }
    int ownerRank() const { return 0; }
    int topologicalDimension() const { return 0; }
    int physicalDimension() const { return 0; }
    void boundingBox( Teuchos::Tuple<double, 6> &bounds ) const {}
    bool inBlock( const int block_id ) const { return false; }
    bool onBoundary( const int boundary_id ) const { return false; }
    Teuchos::RCP<DataTransferKit::EntityExtraData> extraData() const
    {
        return Teuchos::null;
    }

  private:
    DataTransferKit::EntityId d_id;
};

//---------------------------------------------------------------------------//
// Entity Implementation.
//---------------------------------------------------------------------------//
class TestEntity : public DataTransferKit::Entity
{
  public:
    TestEntity( int id )
    {
        this->b_entity_impl = Teuchos::rcp( new TestEntityImpl( id ) );
    }
};

//---------------------------------------------------------------------------//
// Helper predicates.
//---------------------------------------------------------------------------//
std::function<bool( DataTransferKit::Entity )> even_func =
    []( DataTransferKit::Entity e ) { return ( ( e.id() % 2 ) == 0 ); };
std::function<bool( DataTransferKit::Entity )> odd_func =
    []( DataTransferKit::Entity e ) { return ( ( e.id() % 2 ) == 1 ); };
std::function<bool( DataTransferKit::Entity )> two_func =
    []( DataTransferKit::Entity e ) { return ( ( e.id() % 10 ) == 2 ); };

//---------------------------------------------------------------------------//
// EntityIterator implementation.
//---------------------------------------------------------------------------//
class VectorIterator : public DataTransferKit::EntityIterator
{
  public:
    /*!
     * \brief Default constructor.
     */
    VectorIterator()
        : d_it_entity( NULL )
    { /* ... */
    }

    /*!
     * \brief Constructor.
     */
    VectorIterator(
        const Teuchos::RCP<std::vector<DataTransferKit::Entity>> &entities )
        : d_entities( entities )
        , d_vec_it( d_entities->begin() )
    {
        if ( d_vec_it != d_entities->end() )
        {
            d_it_entity = &( *d_vec_it );
        }
    }

    /*!
     * \brief Predicate constructor.
     */
    VectorIterator(
        const Teuchos::RCP<std::vector<DataTransferKit::Entity>> &entities,
        const std::function<bool( DataTransferKit::Entity )> &predicate )
        : d_entities( entities )
        , d_vec_it( d_entities->begin() )
    {
        this->b_predicate = predicate;
        if ( d_vec_it != d_entities->end() )
        {
            d_it_entity = &( *d_vec_it );
        }
    }

    /*!
     * \brief Copy constructor.
     */
    VectorIterator( const VectorIterator &rhs )
        : d_entities( rhs.d_entities )
        , d_vec_it( d_entities->begin() +
                    std::distance( rhs.d_entities->begin(), rhs.d_vec_it ) )
    {
        this->b_predicate = rhs.b_predicate;
        if ( d_vec_it != d_entities->end() )
        {
            d_it_entity = &( *d_vec_it );
        }
    }

    /*!
     * \brief Assignment operator.
     */
    VectorIterator &operator=( const VectorIterator &rhs )
    {
        this->b_predicate = rhs.b_predicate;
        if ( &rhs == this )
        {
            return *this;
        }
        this->d_entities = rhs.d_entities;
        this->d_vec_it = this->d_entities->begin() +
                         std::distance( rhs.d_entities->begin(), rhs.d_vec_it );
        if ( d_vec_it != d_entities->end() )
        {
            d_it_entity = &( *d_vec_it );
        }
        return *this;
    }

    // Pre-increment operator.
    DataTransferKit::EntityIterator &operator++() override
    {
        ++d_vec_it;
        return *this;
    }

    // Dereference operator.
    DataTransferKit::Entity &operator*(void)override
    {
        this->operator->();
        return *d_it_entity;
    }

    // Dereference operator.
    DataTransferKit::Entity *operator->(void)override
    {
        d_it_entity = &( *d_vec_it );
        return d_it_entity;
    }

    // Equal comparison operator.
    bool operator==( const DataTransferKit::EntityIterator &rhs ) const override
    {
        const VectorIterator *rhs_vec =
            static_cast<const VectorIterator *>( &rhs );
        const VectorIterator *rhs_vec_impl =
            static_cast<const VectorIterator *>(
                rhs_vec->b_iterator_impl.get() );
        return ( rhs_vec_impl->d_vec_it == d_vec_it );
    }

    // Not equal comparison operator.
    bool operator!=( const DataTransferKit::EntityIterator &rhs ) const override
    {
        return !( operator==( rhs ) );
    }

    // An iterator assigned to the beginning.
    DataTransferKit::EntityIterator begin() const override
    {
        return VectorIterator( d_entities, this->b_predicate );
    }

    // An iterator assigned to the end.
    DataTransferKit::EntityIterator end() const override
    {
        VectorIterator end_it( d_entities, this->b_predicate );
        end_it.d_vec_it = d_entities->end();
        return end_it;
    }

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    std::unique_ptr<DataTransferKit::EntityIterator> clone() const override
    {
        return std::unique_ptr<DataTransferKit::EntityIterator>(
            new VectorIterator( *this ) );
    }

  private:
    // Vector.
    Teuchos::RCP<std::vector<DataTransferKit::Entity>> d_entities;

    // Iterator to vector.
    typename std::vector<DataTransferKit::Entity>::iterator d_vec_it;

    // Pointer to the current entity.
    DataTransferKit::Entity *d_it_entity;
};

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Constructor tests.
TEUCHOS_UNIT_TEST( EntityIterator, empty_iterator_test )
{
    using namespace DataTransferKit;

    EntityIterator empty_iterator;
    EntityIterator empty_begin = empty_iterator.begin();
    EntityIterator empty_end = empty_iterator.end();

    TEST_EQUALITY( 0, empty_iterator.size() );
    TEST_ASSERT( empty_begin == empty_end );

    int pass_counter = 0;
    for ( empty_iterator = empty_begin; empty_iterator != empty_end;
          ++empty_iterator )
    {
        ++pass_counter;
    }

    TEST_EQUALITY( 0, pass_counter );
}

//---------------------------------------------------------------------------//
// Constructor tests.
TEUCHOS_UNIT_TEST( EntityIterator, constructor_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<Entity>> data =
        Teuchos::rcp( new std::vector<Entity>( num_data ) );
    for ( int i = 0; i < num_data; ++i )
    {
        ( *data )[i] = TestEntity( i );
    }

    // Create an iterator over the vector.
    EntityIterator entity_it = VectorIterator( data );

    // Call the copy constructor.
    EntityIterator it_1( entity_it );
    TEST_ASSERT( entity_it == it_1 );
    TEST_ASSERT( entity_it.size() == it_1.size() );
    TEST_ASSERT( entity_it.begin() == it_1.begin() );
    TEST_ASSERT( entity_it.end() == it_1.end() );
    TEST_EQUALITY( entity_it->id(), it_1->id() );
    ++entity_it;
    ++it_1;
    TEST_EQUALITY( entity_it->id(), it_1->id() );
    ++it_1;
    TEST_ASSERT( entity_it != it_1 );

    // Call the assignment operator.
    it_1 = entity_it;
    TEST_ASSERT( entity_it == it_1 );
    TEST_ASSERT( entity_it.size() == it_1.size() );
    TEST_ASSERT( entity_it.begin() == it_1.begin() );
    TEST_ASSERT( entity_it.end() == it_1.end() );
    TEST_EQUALITY( entity_it->id(), it_1->id() );
    ++entity_it;
    ++it_1;
    TEST_EQUALITY( entity_it->id(), it_1->id() );
    ++it_1;
    TEST_ASSERT( entity_it != it_1 );
}

//---------------------------------------------------------------------------//
// Predicate constructor tests.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_constructor_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<Entity>> data =
        Teuchos::rcp( new std::vector<Entity>( num_data ) );
    for ( int i = 0; i < num_data; ++i )
    {
        ( *data )[i] = TestEntity( i );
    }

    // Create an iterator over the odd components vector.
    EntityIterator odd_it = VectorIterator( data, odd_func );

    // Call the copy constructor.
    EntityIterator it_1( odd_it );
    TEST_ASSERT( odd_it == it_1 );
    TEST_ASSERT( odd_it.size() == it_1.size() );
    TEST_ASSERT( odd_it.begin() == it_1.begin() );
    TEST_ASSERT( odd_it.end() == it_1.end() );
    TEST_EQUALITY( odd_it->id(), it_1->id() );
    ++odd_it;
    ++it_1;
    TEST_EQUALITY( odd_it->id(), it_1->id() );
    ++it_1;
    TEST_ASSERT( odd_it != it_1 );

    // Call the assignment operator.
    it_1 = odd_it;
    TEST_ASSERT( odd_it == it_1 );
    TEST_ASSERT( odd_it.size() == it_1.size() );
    TEST_ASSERT( odd_it.begin() == it_1.begin() );
    TEST_ASSERT( odd_it.end() == it_1.end() );
    TEST_EQUALITY( odd_it->id(), it_1->id() );
    ++odd_it;
    ++it_1;
    TEST_EQUALITY( odd_it->id(), it_1->id() );
    ++it_1;
    TEST_ASSERT( odd_it != it_1 );
    TEST_EQUALITY( odd_it.size(), 5 );
    TEST_EQUALITY( it_1.size(), 5 );

    // Create an iterator over the components vector that end in 2.
    EntityIterator two_it = VectorIterator( data, two_func );
    TEST_EQUALITY( two_it.size(), 1 );
    EntityIterator begin_it = two_it.begin();
    EntityIterator end_it = two_it.end();
    for ( two_it = begin_it; two_it != end_it; ++two_it )
    {
        TEST_ASSERT( two_func( *two_it ) );
    }

    // Assign the two iterator to the odd iterator and check that the
    // predicate was passed.
    two_it = odd_it;
    TEST_EQUALITY( two_it.size(), 5 );
    begin_it = two_it.begin();
    end_it = two_it.end();
    for ( two_it = begin_it; two_it != end_it; ++two_it )
    {
        TEST_ASSERT( odd_func( *two_it ) );
    }
}

//---------------------------------------------------------------------------//
// Basic iterator test.
TEUCHOS_UNIT_TEST( EntityIterator, iterator_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<Entity>> data =
        Teuchos::rcp( new std::vector<Entity>( num_data ) );
    for ( int i = 0; i < num_data; ++i )
    {
        ( *data )[i] = TestEntity( std::rand() );
    }

    // Create an iterator over the vector.
    EntityIterator entity_it = VectorIterator( data );

    // Check size.
    TEST_EQUALITY( num_data, Teuchos::as<int>( entity_it.size() ) );

    // Check the beginning and end.
    EntityIterator begin_it = entity_it.begin();
    EntityIterator end_it = entity_it.end();
    TEST_EQUALITY( num_data, Teuchos::as<int>( begin_it.size() ) );
    TEST_EQUALITY( num_data, Teuchos::as<int>( end_it.size() ) );

    // Check the dereference operators.
    TEST_EQUALITY( entity_it->id(), ( *data )[0].id() );
    TEST_EQUALITY( begin_it->id(), ( *data )[0].id() );

    // Check the comparison operators.
    TEST_ASSERT( begin_it == entity_it );
    TEST_ASSERT( end_it != entity_it );

    // Check the iterator in a for loop.
    std::vector<Entity>::const_iterator data_it;
    for ( entity_it = begin_it, data_it = data->begin(); entity_it != end_it;
          ++entity_it, ++data_it )
    {
        TEST_EQUALITY( data_it->id(), entity_it->id() );
    }

    // Check the increment operators.
    entity_it = begin_it;
    EntityIterator cp_it = entity_it;
    ++entity_it;
    TEST_EQUALITY( entity_it->id(), ( *data )[1].id() );
    TEST_EQUALITY( cp_it->id(), ( *data )[0].id() );
    entity_it++;
    TEST_EQUALITY( entity_it->id(), ( *data )[2].id() );
    TEST_EQUALITY( cp_it->id(), ( *data )[0].id() );
    DataTransferKit::Entity entity = *( ++entity_it );
    TEST_EQUALITY( entity.id(), ( *data )[3].id() );
    TEST_EQUALITY( cp_it->id(), ( *data )[0].id() );
    entity = *( entity_it++ );
    TEST_EQUALITY( entity.id(), ( *data )[3].id() );
    TEST_EQUALITY( entity_it->id(), ( *data )[4].id() );
    TEST_EQUALITY( cp_it->id(), ( *data )[0].id() );
}

//---------------------------------------------------------------------------//
// Single predicate test.
TEUCHOS_UNIT_TEST( EntityIterator, single_predicate_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<Entity>> data =
        Teuchos::rcp( new std::vector<Entity>( num_data ) );
    for ( int i = 0; i < num_data; ++i )
    {
        ( *data )[i] = TestEntity( i );
    }

    // Create an iterator over the vector for the even numbers.
    EntityIterator even_it = VectorIterator( data, even_func );
    TEST_EQUALITY( 5, Teuchos::as<int>( even_it.size() ) );

    // Check the iterator in a for loop.
    EntityIterator begin_it = even_it.begin();
    EntityIterator end_it = even_it.end();
    for ( even_it = begin_it; even_it != end_it; ++even_it )
    {
        TEST_ASSERT( even_func( *even_it ) );
    }

    // Create an iterator over the vector for the odd numbers.
    EntityIterator odd_it = VectorIterator( data, odd_func );
    TEST_EQUALITY( 5, Teuchos::as<int>( odd_it.size() ) );

    // Check the iterator in a for loop.
    begin_it = odd_it.begin();
    end_it = odd_it.end();
    for ( odd_it = begin_it; odd_it != end_it; ++odd_it )
    {
        TEST_ASSERT( odd_func( *odd_it ) );
    }

    // Create an iterator over the vector for numbers with 2 as the last
    // number.
    EntityIterator two_it = VectorIterator( data, two_func );
    TEST_EQUALITY( two_it.size(), 1 );
    begin_it = two_it.begin();
    end_it = two_it.end();
    for ( two_it = begin_it; two_it != end_it; ++two_it )
    {
        TEST_ASSERT( two_func( *two_it ) );
    }
}

//---------------------------------------------------------------------------//
// Predicate intersection test.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_intersection_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<Entity>> data =
        Teuchos::rcp( new std::vector<Entity>( num_data ) );
    for ( int i = 0; i < num_data; ++i )
    {
        ( *data )[i] = TestEntity( i );
    }

    // Create an iterator over the vector for the even and numbers.
    EntityIterator even_odd_intersect_it = VectorIterator(
        data, PredicateComposition::And( even_func, odd_func ) );
    TEST_EQUALITY( even_odd_intersect_it.size(), 0 );
    EntityIterator odd_even_intersect_it = VectorIterator(
        data, PredicateComposition::And( odd_func, even_func ) );
    TEST_EQUALITY( odd_even_intersect_it.size(), 0 );

    // Create the intersection of the even and two set.
    EntityIterator even_two_intersect_it = VectorIterator(
        data, PredicateComposition::And( even_func, two_func ) );
    TEST_EQUALITY( even_two_intersect_it.size(), 1 );
    EntityIterator two_even_intersect_it = VectorIterator(
        data, PredicateComposition::And( two_func, even_func ) );
    TEST_EQUALITY( two_even_intersect_it.size(), 1 );

    // Create the intersection of the two and odd set.
    EntityIterator two_odd_intersection_it =
        VectorIterator( data, PredicateComposition::And( two_func, odd_func ) );
    TEST_EQUALITY( two_odd_intersection_it.size(), 0 );
    EntityIterator odd_two_intersection_it =
        VectorIterator( data, PredicateComposition::And( odd_func, two_func ) );
    TEST_EQUALITY( odd_two_intersection_it.size(), 0 );

    // Intersect the odd set with itself.
    EntityIterator odd_odd_intersect_it =
        VectorIterator( data, PredicateComposition::And( odd_func, odd_func ) );
    TEST_EQUALITY( odd_odd_intersect_it.size(), 5 );
}

//---------------------------------------------------------------------------//
// Predicate union test.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_union_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<Entity>> data =
        Teuchos::rcp( new std::vector<Entity>( num_data ) );
    for ( int i = 0; i < num_data; ++i )
    {
        ( *data )[i] = TestEntity( i );
    }

    // Create the union of the even and odd set.
    EntityIterator even_odd_union_it =
        VectorIterator( data, PredicateComposition::Or( even_func, odd_func ) );
    TEST_EQUALITY( even_odd_union_it.size(), 10 );
    EntityIterator odd_even_union_it =
        VectorIterator( data, PredicateComposition::Or( odd_func, even_func ) );
    TEST_EQUALITY( odd_even_union_it.size(), 10 );

    // Create the union of the even and two set.
    EntityIterator even_two_union_it =
        VectorIterator( data, PredicateComposition::Or( even_func, two_func ) );
    TEST_EQUALITY( even_two_union_it.size(), 5 );
    EntityIterator two_even_union_it =
        VectorIterator( data, PredicateComposition::Or( two_func, even_func ) );
    TEST_EQUALITY( two_even_union_it.size(), 5 );

    // Union the odd set with itself.
    EntityIterator odd_odd_union_it =
        VectorIterator( data, PredicateComposition::Or( odd_func, odd_func ) );
    TEST_EQUALITY( odd_odd_union_it.size(), 5 );

    // Create the union of the two and odd set.
    EntityIterator two_odd_union_it =
        VectorIterator( data, PredicateComposition::Or( two_func, odd_func ) );
    TEST_EQUALITY( two_odd_union_it.size(), 6 );
    EntityIterator odd_two_union_it =
        VectorIterator( data, PredicateComposition::Or( odd_func, two_func ) );
    TEST_EQUALITY( odd_two_union_it.size(), 6 );
}

//---------------------------------------------------------------------------//
// Predicate subtraction test.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_subtraction_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<Entity>> data =
        Teuchos::rcp( new std::vector<Entity>( num_data ) );
    for ( int i = 0; i < num_data; ++i )
    {
        ( *data )[i] = TestEntity( i );
    }

    // Create the subtraction of the even and odd set.
    EntityIterator even_odd_subtraction_it = VectorIterator(
        data, PredicateComposition::AndNot( even_func, odd_func ) );
    TEST_EQUALITY( even_odd_subtraction_it.size(), 5 );
    EntityIterator odd_even_subtraction_it = VectorIterator(
        data, PredicateComposition::AndNot( odd_func, even_func ) );
    TEST_EQUALITY( odd_even_subtraction_it.size(), 5 );

    // Create the subtraction of the even and two set.
    EntityIterator even_two_subtraction_it = VectorIterator(
        data, PredicateComposition::AndNot( even_func, two_func ) );
    TEST_EQUALITY( even_two_subtraction_it.size(), 4 );
    EntityIterator two_even_subtraction_it = VectorIterator(
        data, PredicateComposition::AndNot( two_func, even_func ) );
    TEST_EQUALITY( two_even_subtraction_it.size(), 0 );

    // Subtraction the odd set with itself.
    EntityIterator odd_odd_subtraction_it = VectorIterator(
        data, PredicateComposition::AndNot( odd_func, odd_func ) );
    TEST_EQUALITY( odd_odd_subtraction_it.size(), 0 );

    // Create the subtraction of the odd and two set.
    EntityIterator odd_two_subtraction_it = VectorIterator(
        data, PredicateComposition::AndNot( odd_func, two_func ) );
    TEST_EQUALITY( odd_two_subtraction_it.size(), 5 );
    EntityIterator two_odd_subtraction_it = VectorIterator(
        data, PredicateComposition::AndNot( two_func, odd_func ) );
    TEST_EQUALITY( two_odd_subtraction_it.size(), 1 );
}

//---------------------------------------------------------------------------//
// end tstEntityIterator.cpp
//---------------------------------------------------------------------------//
