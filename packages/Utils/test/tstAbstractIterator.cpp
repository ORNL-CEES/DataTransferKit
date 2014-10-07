//---------------------------------------------------------------------------//
/*!
 * \file tstAbstractIterator.cpp
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

#include <DTK_AbstractIterator.hpp>

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
// AbstractIterator implementation.
//---------------------------------------------------------------------------//
template<class T>
class VectorIterator : public DataTransferKit::AbstractIterator<T>
{
  public:

    /*!
     * \brief Default constructor.
     */
    VectorIterator()
	: d_it_value( 0 )
    {
	this->b_iterator_impl = NULL;
    }

    /*!
     * \brief Constructor.
     */
    VectorIterator( const Teuchos::RCP<std::vector<T> >& values ) 
	: d_values( values )
	, d_vec_it( d_values->begin() )
	, d_it_value( 0 )
    {
	this->b_iterator_impl = NULL;
    }

    /*!
     * \brief Predicate constructor.
     */
    VectorIterator( const Teuchos::RCP<std::vector<T> >& values,
		    const std::function<bool(T)>& predicate ) 
	: d_values( values )
	, d_vec_it( d_values->begin() )
	, d_it_value( 0 )
    {
	this->b_iterator_impl = NULL;
	this->b_predicate = predicate;
    }

    /*!
     * \brief Copy constructor.
     */
    VectorIterator( const VectorIterator<T>& rhs )
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
    VectorIterator& operator=( const VectorIterator<T>& rhs )
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
    DataTransferKit::AbstractIterator<T>& operator++()
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
    bool operator==( const DataTransferKit::AbstractIterator<T>& rhs ) const
    { 
	const VectorIterator<T>* rhs_vec = 
	    static_cast<const VectorIterator<T>*>(&rhs);
	const VectorIterator<T>* rhs_vec_impl = 
	    static_cast<const VectorIterator<T>*>(rhs_vec->b_iterator_impl);
	return ( rhs_vec_impl->d_vec_it == d_vec_it );
    }

    // Not equal comparison operator.
    bool operator!=( const DataTransferKit::AbstractIterator<T>& rhs ) const
    {
	const VectorIterator<T>* rhs_vec = 
	    static_cast<const VectorIterator<T>*>(&rhs);
	const VectorIterator<T>* rhs_vec_impl = 
	    static_cast<const VectorIterator<T>*>(rhs_vec->b_iterator_impl);
	return ( rhs_vec_impl->d_vec_it != d_vec_it );
    }

    // Size of the iterator.
    std::size_t size() const
    { 
	return d_values->size();
    }

    // An iterator assigned to the beginning.
    DataTransferKit::AbstractIterator<T> begin() const
    { 
	return VectorIterator( d_values, this->b_predicate );
    }

    // An iterator assigned to the end.
    DataTransferKit::AbstractIterator<T> end() const
    {
	VectorIterator end_it( d_values, this->b_predicate );
	end_it.d_vec_it = d_values->end();
	return end_it;
    }

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    DataTransferKit::AbstractIterator<T>* clone() const
    {
	return new VectorIterator(*this);
    }

  private:

    // Vector.
    Teuchos::RCP<std::vector<T> > d_values;

    // Iterator To vector.
    typename std::vector<T>::iterator d_vec_it;

    // Pointer to the current value.
    T* d_it_value;
};

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Basic iterator test.
TEUCHOS_UNIT_TEST( AbstractIterator, iterator_test )
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
    AbstractIterator<int> abstract_it = VectorIterator<int>( data );

    // Check size.
    TEST_EQUALITY( num_data, Teuchos::as<int>(abstract_it.size()) );

    // Check the beginning and end.
    AbstractIterator<int> begin_it = abstract_it.begin();
    AbstractIterator<int> end_it = abstract_it.end();
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
    ++abstract_it;
    TEST_EQUALITY( *abstract_it, (*data)[1] );
    abstract_it++;
    TEST_EQUALITY( *abstract_it, (*data)[2] );
    int value = *(++abstract_it);
    TEST_EQUALITY( value, (*data)[3] );
    value = (*abstract_it++);
    TEST_EQUALITY( value, (*data)[3] );
    TEST_EQUALITY( *abstract_it, (*data)[4] );
}

//---------------------------------------------------------------------------//
// Even predicate test.
TEUCHOS_UNIT_TEST( AbstractIterator, even_predicate_test )
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
    AbstractIterator<int> abstract_it = VectorIterator<int>( data );
}
//---------------------------------------------------------------------------//
// end tstAbstractIterator.cpp
//---------------------------------------------------------------------------//

