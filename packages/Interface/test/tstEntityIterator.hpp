//---------------------------------------------------------------------------//
/*!
 * \file tstEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief Abstract iterator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <functional>

#include <DTK_EntityIterator.hpp>
#include <DTK_PredicateComposition.hpp>

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
// Helper predicates.
//---------------------------------------------------------------------------//
std::function<bool(int&)> even_func = [](int& n){ return ((n%2) == 0); };
std::function<bool(int&)> odd_func = [](int& n){ return ((n%2) == 1); };
std::function<bool(int&)> two_func = [](int& n){ return ((n%10) == 2); };

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
	: d_it_value( NULL )
    {
	this->b_iterator_impl = NULL;
    }

    /*!
     * \brief Constructor.
     */
    VectorIterator( const Teuchos::RCP<std::vector<DataTransferKit::Entity> >& values ) 
	: d_values( values )
	, d_vec_it( d_values->begin() )
	, d_it_value( &(*d_vec_it) )
    {
	this->b_iterator_impl = NULL;
    }

    /*!
     * \brief Predicate constructor.
     */
    VectorIterator( const Teuchos::RCP<std::vector<DataTransferKit::Entity> >& values,
		    const std::function<bool(T&)>& predicate ) 
	: d_values( values )
	, d_vec_it( d_values->begin() )
	, d_it_value( &(*d_vec_it) )
    {
	this->b_iterator_impl = NULL;
	this->b_predicate = predicate;
    }

    /*!
     * \brief Copy constructor.
     */
    VectorIterator( const VectorIterator<DataTransferKit::Entity>& rhs )
	: d_values( rhs.d_values )
	, d_vec_it( d_values->begin() + 
		    std::distance(rhs.d_values->begin(),rhs.d_vec_it) )
	, d_it_value( &(*d_vec_it) )
    {
	this->b_iterator_impl = NULL;
	this->b_predicate = rhs.b_predicate;
    }

    /*!
     * \brief Assignment operator.
     */
    VectorIterator& operator=( const VectorIterator<DataTransferKit::Entity>& rhs )
    {
	this->b_iterator_impl = NULL;
	this->b_predicate = rhs.b_predicate;
	if ( &rhs == this )
	{
	    return *this;
	}
	this->d_values = rhs.d_values;
	this->d_vec_it = this->d_values->begin() + 
			 std::distance(rhs.d_values->begin(),rhs.d_vec_it);
	this->d_it_value = &(*(this->d_vec_it));
	return *this;
    }

    /*!
     * \brief Destructor.
     */
    ~VectorIterator()
    {
	this->b_iterator_impl = NULL;
    }

    // Pre-increment operator.
    DataTransferKit::EntityIterator& operator++()
    {
	++d_vec_it;
	return *this;
    }

    // Dereference operator.
    T& operator*(void)
    {
	this->operator->();
	return *d_it_value;
    }

    // Dereference operator.
    T* operator->(void)
    {
	d_it_value = &(*d_vec_it);
	return d_it_value;
    }

    // Equal comparison operator.
    bool operator==( const DataTransferKit::EntityIterator& rhs ) const
    { 
	const VectorIterator<DataTransferKit::Entity>* rhs_vec = 
	    static_cast<const VectorIterator<DataTransferKit::Entity>*>(&rhs);
	const VectorIterator<DataTransferKit::Entity>* rhs_vec_impl = 
	    static_cast<const VectorIterator<DataTransferKit::Entity>*>(rhs_vec->b_iterator_impl);
	return ( rhs_vec_impl->d_vec_it == d_vec_it );
    }

    // Not equal comparison operator.
    bool operator!=( const DataTransferKit::EntityIterator& rhs ) const
    {
	return !( operator==(rhs) );
    }

    // Size of the iterator.
    std::size_t size() const
    { 
	return d_values->size();
    }

    // An iterator assigned to the beginning.
    DataTransferKit::EntityIterator begin() const
    { 
	return VectorIterator( d_values , this->b_predicate );
    }

    // An iterator assigned to the end.
    DataTransferKit::EntityIterator end() const
    {
	VectorIterator end_it( d_values, this->b_predicate );
	end_it.d_vec_it = d_values->end();
	return end_it;
    }

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    DataTransferKit::EntityIterator* clone() const
    {
	return new VectorIterator(*this);
    }

  private:

    // Vector.
    Teuchos::RCP<std::vector<DataTransferKit::Entity> > d_values;

    // Iterator To vector.
    typename std::vector<DataTransferKit::Entity>::iterator d_vec_it;

    // Pointer to the current entity.
    Entity* d_it_value;
};

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Constructor tests.
TEUCHOS_UNIT_TEST( EntityIterator, constructor_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<int> > data =
	Teuchos::rcp( new std::vector<int>(num_data) );
    for ( int i = 0; i < num_data; ++i )
    {
	(*data)[i] = i;
    }

    // Create an iterator over the vector.
    EntityIterator<int> abstract_it = VectorIterator<int>( data );

    // Call the copy constructor.
    EntityIterator<int> it_1( abstract_it );
    TEST_ASSERT( abstract_it == it_1 );
    TEST_ASSERT( abstract_it.size() == it_1.size() );
    TEST_ASSERT( abstract_it.begin() == it_1.begin() );
    TEST_ASSERT( abstract_it.end() == it_1.end() );
    TEST_EQUALITY( *abstract_it, *it_1 );
    ++abstract_it;
    ++it_1;
    TEST_EQUALITY( *abstract_it, *it_1 );
    ++it_1;
    TEST_ASSERT( abstract_it != it_1 );

    // Call the assignment operator.
    it_1 = abstract_it;
    TEST_ASSERT( abstract_it == it_1 );
    TEST_ASSERT( abstract_it.size() == it_1.size() );
    TEST_ASSERT( abstract_it.begin() == it_1.begin() );
    TEST_ASSERT( abstract_it.end() == it_1.end() );
    TEST_EQUALITY( *abstract_it, *it_1 );
    ++abstract_it;
    ++it_1;
    TEST_EQUALITY( *abstract_it, *it_1 );
    ++it_1;
    TEST_ASSERT( abstract_it != it_1 );
}

//---------------------------------------------------------------------------//
// Predicate constructor tests.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_constructor_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<int> > data =
	Teuchos::rcp( new std::vector<int>(num_data) );
    for ( int i = 0; i < num_data; ++i )
    {
	(*data)[i] = i;
    }

    // Create an iterator over the odd components vector.
    EntityIterator<int> odd_it = VectorIterator<int>( data, odd_func );

    // Call the copy constructor.
    EntityIterator<int> it_1( odd_it );
    TEST_ASSERT( odd_it == it_1 );
    TEST_ASSERT( odd_it.size() == it_1.size() );
    TEST_ASSERT( odd_it.begin() == it_1.begin() );
    TEST_ASSERT( odd_it.end() == it_1.end() );
    TEST_EQUALITY( *odd_it, *it_1 );
    ++odd_it;
    ++it_1;
    TEST_EQUALITY( *odd_it, *it_1 );
    ++it_1;
    TEST_ASSERT( odd_it != it_1 );

    // Call the assignment operator.
    it_1 = odd_it;
    TEST_ASSERT( odd_it == it_1 );
    TEST_ASSERT( odd_it.size() == it_1.size() );
    TEST_ASSERT( odd_it.begin() == it_1.begin() );
    TEST_ASSERT( odd_it.end() == it_1.end() );
    TEST_EQUALITY( *odd_it, *it_1 );
    ++odd_it;
    ++it_1;
    TEST_EQUALITY( *odd_it, *it_1 );
    ++it_1;
    TEST_ASSERT( odd_it != it_1 );
    TEST_EQUALITY( odd_it.size(), 5 );
    TEST_EQUALITY( it_1.size(), 5 );

    // Create an iterator over the components vector that end in 2.
    EntityIterator<int> two_it = VectorIterator<int>( data, two_func );
    TEST_EQUALITY( two_it.size(), 1 );
    EntityIterator<int>begin_it = two_it.begin();
    EntityIterator<int> end_it = two_it.end();
    for ( two_it = begin_it; two_it != end_it; ++two_it )
    {
    	TEST_ASSERT( two_func(*two_it) );
    }

    // Assign the two iterator to the odd iterator and check that the
    // predicate was passed.
    two_it = odd_it;
    TEST_EQUALITY( two_it.size(), 5 );
    begin_it = two_it.begin();
    end_it = two_it.end();
    for ( two_it = begin_it; two_it != end_it; ++two_it )
    {
    	TEST_ASSERT( odd_func(*two_it) );
    }
}

//---------------------------------------------------------------------------//
// Basic iterator test.
TEUCHOS_UNIT_TEST( EntityIterator, iterator_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<int> > data =
	Teuchos::rcp( new std::vector<int>(num_data) );
    for ( int i = 0; i < num_data; ++i )
    {
	(*data)[i] = std::rand();
    }

    // Create an iterator over the vector.
    EntityIterator<int> abstract_it = VectorIterator<int>( data );

    // Check size.
    TEST_EQUALITY( num_data, Teuchos::as<int>(abstract_it.size()) );

    // Check the beginning and end.
    EntityIterator<int> begin_it = abstract_it.begin();
    EntityIterator<int> end_it = abstract_it.end();
    TEST_EQUALITY( num_data, Teuchos::as<int>(begin_it.size()) );
    TEST_EQUALITY( num_data, Teuchos::as<int>(end_it.size()) );

    // Check the dereference operators.
    TEST_EQUALITY( *abstract_it, (*data)[0] );
    TEST_EQUALITY( *begin_it, (*data)[0] );

    // Check the comparison operators.
    TEST_ASSERT( begin_it == abstract_it );
    TEST_ASSERT( end_it != abstract_it );

    // Check the iterator in a for loop.
    std::vector<int>::const_iterator data_it;
    for ( abstract_it = begin_it, data_it = data->begin(); 
    	  abstract_it != end_it; 
    	  ++abstract_it, ++data_it )
    {
    	TEST_EQUALITY( *data_it, *abstract_it );
    }

    // Check the increment operators.
    abstract_it = begin_it;
    EntityIterator<int> cp_it = abstract_it;
    ++abstract_it;
    TEST_EQUALITY( *abstract_it, (*data)[1] );
    TEST_EQUALITY( *cp_it, (*data)[0] );
    abstract_it++;
    TEST_EQUALITY( *abstract_it, (*data)[2] );
    TEST_EQUALITY( *cp_it, (*data)[0] );
    int value = *(++abstract_it);
    TEST_EQUALITY( value, (*data)[3] );
    TEST_EQUALITY( *cp_it, (*data)[0] );
    value = (*abstract_it++);
    TEST_EQUALITY( value, (*data)[3] );
    TEST_EQUALITY( *abstract_it, (*data)[4] );
    TEST_EQUALITY( *cp_it, (*data)[0] );
}

//---------------------------------------------------------------------------//
// Single predicate test.
TEUCHOS_UNIT_TEST( EntityIterator, single_predicate_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<int> > data =
    	Teuchos::rcp( new std::vector<int>(num_data) );
    for ( int i = 0; i < num_data; ++i )
    {
    	(*data)[i] = i;
    }

    // Create an iterator over the vector for the even numbers.
    EntityIterator<int> even_it = VectorIterator<int>( data, even_func );
    TEST_EQUALITY( 5, Teuchos::as<int>(even_it.size()) );

    // Check the iterator in a for loop.
    EntityIterator<int> begin_it = even_it.begin();
    EntityIterator<int> end_it = even_it.end();
    for ( even_it = begin_it;
    	  even_it != end_it; 
    	  ++even_it )
    {
    	TEST_ASSERT( even_func(*even_it) );
    }

    // Create an iterator over the vector for the odd numbers.
    EntityIterator<int> odd_it = VectorIterator<int>( data, odd_func );
    TEST_EQUALITY( 5, Teuchos::as<int>(odd_it.size()) );

    // Check the iterator in a for loop.
    begin_it = odd_it.begin();
    end_it = odd_it.end();
    for ( odd_it = begin_it;
    	  odd_it != end_it; 
    	  ++odd_it )
    {
    	TEST_ASSERT( odd_func(*odd_it) );
    }

    // Create an iterator over the vector for numbers with 2 as the last
    // number.
    EntityIterator<int> two_it = VectorIterator<int>( data, two_func );
    TEST_EQUALITY( two_it.size(), 1 );
    begin_it = two_it.begin();
    end_it = two_it.end();
    for ( two_it = begin_it; two_it != end_it; ++two_it )
    {
    	TEST_ASSERT( two_func(*two_it) );
    }
}

//---------------------------------------------------------------------------//
// Predicate intersection test.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_intersection_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<int> > data =
    	Teuchos::rcp( new std::vector<int>(num_data) );
    for ( int i = 0; i < num_data; ++i )
    {
    	(*data)[i] = i;
    }

    // Create an iterator over the vector for the even and numbers.
    EntityIterator<int> even_odd_intersect_it = 
	VectorIterator<int>( data, PredicateComposition::And(even_func,odd_func) );
    TEST_EQUALITY( even_odd_intersect_it.size(), 0 );
    EntityIterator<int> odd_even_intersect_it = 
	VectorIterator<int>( data, PredicateComposition::And(odd_func,even_func) );
    TEST_EQUALITY( odd_even_intersect_it.size(), 0 );

    // Create the intersection of the even and two set.
    EntityIterator<int> even_two_intersect_it = 
	VectorIterator<int>( data, PredicateComposition::And(even_func,two_func) );
    TEST_EQUALITY( even_two_intersect_it.size(), 1 );
    EntityIterator<int> two_even_intersect_it = 
	VectorIterator<int>( data, PredicateComposition::And(two_func,even_func) );
    TEST_EQUALITY( two_even_intersect_it.size(), 1 );

    // Create the intersection of the two and odd set.
    EntityIterator<int> two_odd_intersection_it = 
	VectorIterator<int>( data, PredicateComposition::And(two_func,odd_func) );
    TEST_EQUALITY( two_odd_intersection_it.size(), 0 );
    EntityIterator<int> odd_two_intersection_it = 
	VectorIterator<int>( data, PredicateComposition::And(odd_func,two_func) );
    TEST_EQUALITY( odd_two_intersection_it.size(), 0 );

    // Intersect the odd set with itself.
    EntityIterator<int> odd_odd_intersect_it = 
	VectorIterator<int>( data, PredicateComposition::And(odd_func,odd_func) );
    TEST_EQUALITY( odd_odd_intersect_it.size(), 5 );
}

//---------------------------------------------------------------------------//
// Predicate union test.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_union_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<int> > data =
    	Teuchos::rcp( new std::vector<int>(num_data) );
    for ( int i = 0; i < num_data; ++i )
    {
    	(*data)[i] = i;
    }

    // Create the union of the even and odd set.
    EntityIterator<int> even_odd_union_it = 
	VectorIterator<int>( data, PredicateComposition::Or(even_func,odd_func) );
    TEST_EQUALITY( even_odd_union_it.size(), 10 );
    EntityIterator<int> odd_even_union_it = 
	VectorIterator<int>( data, PredicateComposition::Or(odd_func,even_func) );
    TEST_EQUALITY( odd_even_union_it.size(), 10 );

    // Create the union of the even and two set.
    EntityIterator<int> even_two_union_it = 
	VectorIterator<int>( data, PredicateComposition::Or(even_func,two_func) );
    TEST_EQUALITY( even_two_union_it.size(), 5 );
    EntityIterator<int> two_even_union_it = 
	VectorIterator<int>( data, PredicateComposition::Or(two_func,even_func) );
    TEST_EQUALITY( two_even_union_it.size(), 5 );

    // Union the odd set with itself.
    EntityIterator<int> odd_odd_union_it = 
	VectorIterator<int>( data, PredicateComposition::Or(odd_func,odd_func) );
    TEST_EQUALITY( odd_odd_union_it.size(), 5 );

    // Create the union of the two and odd set.
    EntityIterator<int> two_odd_union_it = 
	VectorIterator<int>( data, PredicateComposition::Or(two_func,odd_func) );
    TEST_EQUALITY( two_odd_union_it.size(), 6 );
    EntityIterator<int> odd_two_union_it = 
	VectorIterator<int>( data, PredicateComposition::Or(odd_func,two_func) );
    TEST_EQUALITY( odd_two_union_it.size(), 6 );
}

//---------------------------------------------------------------------------//
// Predicate subtraction test.
TEUCHOS_UNIT_TEST( EntityIterator, predicate_subtraction_test )
{
    using namespace DataTransferKit;

    // Create a vector.
    int num_data = 10;
    Teuchos::RCP<std::vector<int> > data =
    	Teuchos::rcp( new std::vector<int>(num_data) );
    for ( int i = 0; i < num_data; ++i )
    {
    	(*data)[i] = i;
    }

    // Create the subtraction of the even and odd set.
    EntityIterator<int> even_odd_subtraction_it = 
	VectorIterator<int>( data, PredicateComposition::AndNot(even_func,odd_func) );
    TEST_EQUALITY( even_odd_subtraction_it.size(), 5 );
    EntityIterator<int> odd_even_subtraction_it = 
	VectorIterator<int>( data, PredicateComposition::AndNot(odd_func,even_func) );
    TEST_EQUALITY( odd_even_subtraction_it.size(), 5 );

    // Create the subtraction of the even and two set.
    EntityIterator<int> even_two_subtraction_it = 
	VectorIterator<int>( data, PredicateComposition::AndNot(even_func,two_func) );
    TEST_EQUALITY( even_two_subtraction_it.size(), 4 );
    EntityIterator<int> two_even_subtraction_it = 
	VectorIterator<int>( data, PredicateComposition::AndNot(two_func,even_func) );
    TEST_EQUALITY( two_even_subtraction_it.size(), 0 );

    // Subtraction the odd set with itself.
    EntityIterator<int> odd_odd_subtraction_it = 
	VectorIterator<int>( data, PredicateComposition::AndNot(odd_func,odd_func) );
    TEST_EQUALITY( odd_odd_subtraction_it.size(), 0 );

    // Create the subtraction of the odd and two set.
    EntityIterator<int> odd_two_subtraction_it = 
	VectorIterator<int>( data, PredicateComposition::AndNot(odd_func,two_func) );
    TEST_EQUALITY( odd_two_subtraction_it.size(), 5 );
    EntityIterator<int> two_odd_subtraction_it = 
	VectorIterator<int>( data, PredicateComposition::AndNot(two_func,odd_func) );
    TEST_EQUALITY( two_odd_subtraction_it.size(), 1 );
}

//---------------------------------------------------------------------------//
// end tstEntityIterator.cpp
//---------------------------------------------------------------------------//

