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

#include <DTK_Classic_FieldManager.hpp>
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
// Field implementation.
//---------------------------------------------------------------------------//
class ArrayField
{
  public:

    typedef double value_type;
    typedef Teuchos::Array<double>::size_type size_type;
    typedef Teuchos::Array<double>::iterator iterator;
    typedef Teuchos::Array<double>::const_iterator const_iterator;

    ArrayField( size_type size, int dim )
	: d_dim( dim )
	, d_data( size )
    { /* ... */ }

    ~ArrayField()
    { /* ... */ }

    int dim() const
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
    int d_dim;
    Teuchos::Array<double> d_data;
};

//---------------------------------------------------------------------------//
// DTK implementations.
//---------------------------------------------------------------------------//
namespace DataTransferKit
{
namespace Classic
{
//---------------------------------------------------------------------------//
// Field Traits specification for ArrayField
template<>
class FieldTraits<ArrayField>
{
  public:

    typedef ArrayField                    field_type;
    typedef ArrayField::value_type        value_type;
    typedef ArrayField::size_type         size_type;
    typedef ArrayField::iterator          iterator;
    typedef ArrayField::const_iterator    const_iterator;

    static inline size_type dim( const ArrayField& field )
    { return field.dim(); }

    static inline size_type size( const ArrayField& field )
    { return field.size(); }

    static inline bool empty( const ArrayField& field )
    { return field.empty(); }

    static inline iterator begin( ArrayField& field )
    { return field.begin(); }

    static inline const_iterator begin( const ArrayField& field )
    { return field.begin(); }

    static inline iterator end( ArrayField& field )
    { return field.end(); }

    static inline const_iterator end( const ArrayField& field )
    { return field.end(); }
};

} // end namespace Classic
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Unit tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( FieldManager, field_manager_test )
{
    using namespace DataTransferKit::Classic;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Setup field manager.
    Teuchos::RCP<ArrayField> array_field = 
	Teuchos::rcp( new ArrayField( 12, 4 ) );
    FieldManager<ArrayField> field_manager( array_field, comm );

    // Check the field manager.
    Teuchos::RCP<ArrayField> manager_field = field_manager.field();
    TEST_ASSERT( FieldTraits<ArrayField>::dim( *manager_field ) == 
		 FieldTraits<ArrayField>::dim( *array_field ) );
    TEST_ASSERT( FieldTraits<ArrayField>::size( *manager_field ) == 
		 FieldTraits<ArrayField>::size( *array_field ) );
    TEST_ASSERT( field_manager.comm() == comm );
}

//---------------------------------------------------------------------------//
// end tstFieldManager.cpp
//---------------------------------------------------------------------------//
