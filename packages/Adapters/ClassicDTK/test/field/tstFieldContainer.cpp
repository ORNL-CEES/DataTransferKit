//---------------------------------------------------------------------------//
/*!
 * \file tstFieldManager.cpp
 * \author Stuart R. Slattery
 * \brief FieldManager unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_Classic_FieldContainer.hpp>
#include <DTK_Classic_FieldTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_TypeTraits.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Unit tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( FieldManager, field_manager_test )
{
    using namespace DataTransferKit::Classic;
    typedef FieldTraits< FieldContainer<double> >  FT;

    // Setup some data.
    int field_size = 20;
    Teuchos::ArrayRCP<double> my_data( field_size );
    Teuchos::ArrayRCP<double>::iterator data_iterator;
    double value = 0.0;
    for ( data_iterator = my_data.begin(); data_iterator != my_data.end();
	  ++data_iterator )
    {
	value = 1.0 * std::distance( my_data.begin(), data_iterator );
	*data_iterator = value;
    }

    // Setup field container.
    int field_dimension = 4;
    FieldContainer<double> container( my_data, field_dimension );
    
    // Check the field container.
    TEST_ASSERT( container.dim() == field_dimension );
    TEST_ASSERT( FT::dim( container ) == field_dimension );
    TEST_ASSERT( container.size() == field_size );
    TEST_ASSERT( FT::size( container ) == field_size );
    TEST_ASSERT( !container.empty() );
    TEST_ASSERT( !FT::empty( container ) );

    FieldContainer<double>::iterator container_iterator;
    for ( container_iterator = container.begin();
	  container_iterator != container.end();
	  ++container_iterator )
    {
	value = std::distance( container.begin(), container_iterator );
	TEST_ASSERT( *container_iterator == value );
    }

    FT::iterator traits_iterator;
    for ( traits_iterator = FT::begin( container );
	  traits_iterator != FT::begin( container );
	  ++traits_iterator )
    {
	value = std::distance( FT::begin( container ), traits_iterator );
	TEST_ASSERT( *traits_iterator == value );
    }

    // Build an empty container.
    Teuchos::ArrayRCP<double> empty_data( 0 );
    FieldContainer<double> empty_container( empty_data, 0 );
    TEST_ASSERT( empty_container.empty() );
    TEST_ASSERT( FT::empty( empty_container ) );
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//
