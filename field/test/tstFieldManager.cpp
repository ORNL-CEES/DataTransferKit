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

#include <DTK_FieldManager.hpp>
#include <DTK_FieldTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
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
// Field implementation.
//---------------------------------------------------------------------------//
class MyField
{
  public:

    typedef double value_type;
    typedef Teuchos::Array<double>::size_type size_type;
    typedef Teuchos::Array<double>::iterator iterator;
    typedef Teuchos::Array<double>::const_iterator const_iterator;

    MyField( size_type size, std::size_t dim )
	: d_dim( dim )
	, d_data( size )
    { /* ... */ }

    ~MyField()
    { /* ... */ }

    std::size_t dim() const
    { return d_dim; }

    size_type size() const
    { return d_data.size(); }

    bool empty() const
    { return d_data.empty(); }

    iterator begin()
    { return d_data.begin(); }

    const_iterator begin() const
    { return d_data.begin(); }

    iterator end()
    { return d_data.end(); }

    const_iterator end() const
    { return d_data.end(); }

    const Teuchos::Array<double>& getData() const
    { return d_data; }

  private:
    std::size_t d_dim;
    Teuchos::Array<double> d_data;
};

//---------------------------------------------------------------------------//
// DTK implementations.
//---------------------------------------------------------------------------//
namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Field Traits specification for MyField
template<>
class FieldTraits<MyField>
{
  public:

    typedef MyField                    field_type;
    typedef double                                    value_type;
    typedef MyField::size_type         size_type;
    typedef MyField::iterator          iterator;
    typedef MyField::const_iterator    const_iterator;

    static inline size_type dim( const MyField& field )
    { return field.dim(); }

    static inline size_type size( const MyField& field )
    { return field.size(); }

    static inline bool empty( const MyField& field )
    { return field.empty(); }

    static inline iterator begin( MyField& field )
    { return field.begin(); }

    static inline const_iterator begin( const MyField& field )
    { return field.begin(); }

    static inline iterator end( MyField& field )
    { return field.end(); }

    static inline const_iterator end( const MyField& field )
    { return field.end(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Unit tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( FieldManager, field_manager_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Setup field manager.
    MyField my_field( 12, 4 );
    FieldManager<MyField> field_manager( my_field, comm );

    // Check the field manager.
    MyField manager_field = field_manager.field();
    TEST_ASSERT( FieldTraits<MyField>::dim( manager_field ) == 
		 FieldTraits<MyField>::dim( my_field ) );
    TEST_ASSERT( FieldTraits<MyField>::size( manager_field ) == 
		 FieldTraits<MyField>::size( my_field ) );
    TEST_ASSERT( field_manager.comm() == comm );
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//
