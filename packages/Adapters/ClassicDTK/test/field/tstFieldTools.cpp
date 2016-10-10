//---------------------------------------------------------------------------//
/*!
 * \file tstFieldTools.cpp
 * \author Stuart R. Slattery
 * \brief Field tools unit test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include <DTK_FieldTools.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Ptr.hpp>

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
// Helper Functions
//---------------------------------------------------------------------------//
bool softEquivalence( double a1, double a2, double tol=1.0e-6 )
{
    if ( std::abs( a1 - a2 ) < tol )
    {
        return true;
    }
    else
    {
        return false;
    }
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

    Teuchos::Array<double>& getData()
    { return d_data; }

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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FieldTools, scalar_field_test )
{
    using namespace DataTransferKit;

    // Setup random numbers for the test.
    int num_rand = 4;
    double rand_max = 10.0;
    Teuchos::Array<double> random_numbers( num_rand );
    std::srand( 1 );
    for ( int i = 0; i < num_rand; ++i )
    {
        random_numbers[i] = rand_max * (double) std::rand() / RAND_MAX;
    }

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup a field manager.
    int field_dim = 1;
    typedef FieldTraits<ArrayField> FT;
    Teuchos::RCP<ArrayField> array_field =
        Teuchos::rcp( new ArrayField( field_dim*(my_rank+1), field_dim ) );
    FieldManager<ArrayField> field_manager( array_field, comm );

    // Test the field tools.
    typedef FieldTools<ArrayField> Tools;

    // Dimension iterators.
    FT::iterator dim_begin, dim_end;
    for ( int d = 0; d < field_dim; ++d )
    {
        dim_begin = Tools::dimBegin( *field_manager.field(), d );
        dim_end = Tools::dimEnd( *field_manager.field(), d );
        TEST_ASSERT( std::distance( dim_begin, dim_end ) == my_rank+1 );
    }

    FT::const_iterator const_dim_begin, const_dim_end;
    for ( int d = 0; d < field_dim; ++d )
    {
        const_dim_begin = Tools::dimBegin( *field_manager.field(), d );
        const_dim_end = Tools::dimEnd( *field_manager.field(), d );
        TEST_ASSERT( std::distance( const_dim_begin, const_dim_end )
                     == my_rank+1 );
    }

    // Length.
    int global_size = 0;
    for ( int i = 0; i < my_size; ++i )
    {
        global_size += field_dim*(i+1);
    }
    TEST_ASSERT(
        Tools::globalSize( *field_manager.field(), field_manager.comm() )
        == global_size );

    // Scalar filling.
    Tools::putScalar( *field_manager.field(), random_numbers[0] );

    // Views.
    double local_val;
    Teuchos::ArrayRCP<const double> const_view =
        Tools::view( *field_manager.field() );
    Teuchos::ArrayRCP<const double>::const_iterator const_view_iterator;
    for ( const_view_iterator = const_view.begin();
          const_view_iterator != const_view.end();
          ++const_view_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *const_view_iterator == local_val );
    }

    Teuchos::ArrayRCP<double> view =
        Tools::nonConstView( *field_manager.field() );
    Teuchos::ArrayRCP<double>::const_iterator view_iterator;
    for ( view_iterator = view.begin();
          view_iterator != view.end();
          ++view_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *view_iterator == local_val );
    }

    Teuchos::ArrayRCP<double> copy =
        Tools::copy( *field_manager.field() );
    Teuchos::ArrayRCP<double>::const_iterator copy_iterator;
    for ( copy_iterator = copy.begin();
          copy_iterator != copy.end();
          ++copy_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *copy_iterator == local_val );
    }

    // Multiple scalar filling.
    Teuchos::Array<double> fillers( field_dim );
    for ( int d = 0; d < field_dim; ++d )
    {
        fillers[d] = (d+1)*random_numbers[1];
    }

    Tools::putScalar( *field_manager.field(), fillers() );

    // Dim iterators.
    FT::const_iterator const_dim_iterator;
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( const_dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              const_dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++const_dim_iterator )
        {
            local_val = (d+1) * random_numbers[1];
            TEST_ASSERT( *const_dim_iterator == local_val );
        }
    }

    FT::iterator dim_iterator;
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1];
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Scaling.
    Tools::scale( *field_manager.field(), random_numbers[2] );
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1] * random_numbers[2];
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Multiple scaling.
    Teuchos::Array<double> multipliers( field_dim );
    for ( int d = 0; d < field_dim; ++d )
    {
        multipliers[d] = -(d+1) - random_numbers[3];
    }
    Tools::scale( *field_manager.field(), multipliers() );
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1]
                        * random_numbers[2] * (-(d+1) - random_numbers[3]);
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Norms.
    Teuchos::Array<double> inf_norms;
    Tools::normInf( *field_manager.field(), field_manager.comm(), inf_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = (d+1) * random_numbers[1] * random_numbers[2]
                    * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( inf_norms[d], std::abs( local_val ) ) );
    }

    Teuchos::Array<double> one_norms;
    Tools::norm1( *field_manager.field(), field_manager.comm(), one_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                    (d+1) * random_numbers[1]
                    * random_numbers[2] * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( one_norms[d], std::abs( local_val ) ) );
    }

    Teuchos::Array<double> two_norms;
    Tools::norm2( *field_manager.field(), field_manager.comm(), two_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                   std::pow(
                       std::abs( (d+1) * random_numbers[1] * random_numbers[2]
                                 * (-(d+1) - random_numbers[3]) ), 2.0 );
        local_val = std::pow( local_val, 1.0/2.0 );
        TEST_ASSERT( softEquivalence( two_norms[d], local_val ) );
    }

    Teuchos::Array<double> q_norms;
    int q = std::floor( (10 * std::rand()) / RAND_MAX ) + 1;
    Tools::normQ( *field_manager.field(), field_manager.comm(), q, q_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                    std::pow(
                        std::abs( (d+1) * random_numbers[1] * random_numbers[2]
                                  * (-(d+1) - random_numbers[3]) ), q );
        local_val = std::pow( local_val, 1.0/q );
        TEST_ASSERT( softEquivalence( q_norms[d], local_val ) );
    }

    // Average.
    Teuchos::Array<double> averages;
    Tools::average( *field_manager.field(), field_manager.comm(), averages );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = (d+1) * random_numbers[1] * random_numbers[2]
                    * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( averages[d], local_val ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FieldTools, vector_field_test )
{
    using namespace DataTransferKit;

    // Setup random numbers for the test.
    int num_rand = 4;
    double rand_max = 10.0;
    Teuchos::Array<double> random_numbers( num_rand );
    std::srand( 1 );
    for ( int i = 0; i < num_rand; ++i )
    {
        random_numbers[i] = rand_max * (double) std::rand() / RAND_MAX;
    }

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup a field manager.
    int field_dim = 3;
    typedef FieldTraits<ArrayField> FT;
    Teuchos::RCP<ArrayField> array_field =
        Teuchos::rcp( new ArrayField( field_dim*(my_rank+1), field_dim ) );
    FieldManager<ArrayField> field_manager( array_field, comm );

    // Test the field tools.
    typedef FieldTools<ArrayField> Tools;

    // Dimension iterators.
    FT::iterator dim_begin, dim_end;
    for ( int d = 0; d < field_dim; ++d )
    {
        dim_begin = Tools::dimBegin( *field_manager.field(), d );
        dim_end = Tools::dimEnd( *field_manager.field(), d );
        TEST_ASSERT( std::distance( dim_begin, dim_end ) == my_rank+1 );
    }

    FT::const_iterator const_dim_begin, const_dim_end;
    for ( int d = 0; d < field_dim; ++d )
    {
        const_dim_begin = Tools::dimBegin( *field_manager.field(), d );
        const_dim_end = Tools::dimEnd( *field_manager.field(), d );
        TEST_ASSERT( std::distance( const_dim_begin, const_dim_end )
                     == my_rank+1 );
    }

    // Length.
    int global_size = 0;
    for ( int i = 0; i < my_size; ++i )
    {
        global_size += field_dim*(i+1);
    }
    TEST_ASSERT(
        Tools::globalSize( *field_manager.field(), field_manager.comm() )
        == global_size );

    // Scalar filling.
    Tools::putScalar( *field_manager.field(), random_numbers[0] );

    // Views.
    double local_val;
    Teuchos::ArrayRCP<const double> const_view =
        Tools::view( *field_manager.field() );
    Teuchos::ArrayRCP<const double>::const_iterator const_view_iterator;
    for ( const_view_iterator = const_view.begin();
          const_view_iterator != const_view.end();
          ++const_view_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *const_view_iterator == local_val );
    }

    Teuchos::ArrayRCP<double> view =
        Tools::nonConstView( *field_manager.field() );
    Teuchos::ArrayRCP<double>::iterator view_iterator;
    for ( view_iterator = view.begin();
          view_iterator != view.end();
          ++view_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *view_iterator == local_val );
    }

    Teuchos::ArrayRCP<double> copy =
        Tools::copy( *field_manager.field() );
    Teuchos::ArrayRCP<double>::const_iterator copy_iterator;
    for ( copy_iterator = copy.begin();
          copy_iterator != copy.end();
          ++copy_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *copy_iterator == local_val );
    }

    // Multiple scalar filling.
    Teuchos::Array<double> fillers( field_dim );
    for ( int d = 0; d < field_dim; ++d )
    {
        fillers[d] = (d+1)*random_numbers[1];
    }
    Tools::putScalar( *field_manager.field(), fillers() );

    // Dim iterators.
    FT::const_iterator const_dim_iterator;
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( const_dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              const_dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++const_dim_iterator )
        {
            local_val = (d+1) * random_numbers[1];
            TEST_ASSERT( *const_dim_iterator == local_val );
        }
    }

    FT::iterator dim_iterator;
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1];
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Scaling.
    Tools::scale( *field_manager.field(), random_numbers[2] );
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1] * random_numbers[2];
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Multiple scaling.
    Teuchos::Array<double> multipliers( field_dim );
    for ( int d = 0; d < field_dim; ++d )
    {
        multipliers[d] = -(d+1) - random_numbers[3];
    }
    Tools::scale( *field_manager.field(), multipliers() );
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1]
                        * random_numbers[2] * (-(d+1) - random_numbers[3]);
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Norms.
    Teuchos::Array<double> inf_norms;
    Tools::normInf( *field_manager.field(), field_manager.comm(), inf_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = (d+1) * random_numbers[1] * random_numbers[2]
                    * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( inf_norms[d], std::abs( local_val ) ) );
    }

    Teuchos::Array<double> one_norms;
    Tools::norm1( *field_manager.field(), field_manager.comm(), one_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                   (d+1) * random_numbers[1]
                   * random_numbers[2] * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( one_norms[d], std::abs( local_val ) ) );
    }

    Teuchos::Array<double> two_norms;
    Tools::norm2( *field_manager.field(), field_manager.comm(), two_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                   std::pow(
                       std::abs( (d+1) * random_numbers[1] * random_numbers[2]
                                 * (-(d+1) - random_numbers[3]) ), 2.0 );
        local_val = std::pow( local_val, 1.0/2.0 );
        TEST_ASSERT( softEquivalence( two_norms[d], local_val ) );
    }

    Teuchos::Array<double> q_norms;
    int q = std::floor( (10 * std::rand()) / RAND_MAX ) + 1;
    Tools::normQ( *field_manager.field(), field_manager.comm(), q, q_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                    std::pow(
                        std::abs( (d+1) * random_numbers[1] * random_numbers[2]
                                  * (-(d+1) - random_numbers[3]) ), q );
        local_val = std::pow( local_val, 1.0/q );
        TEST_ASSERT( softEquivalence( q_norms[d], local_val ) );
    }

    // Average.
    Teuchos::Array<double> averages;
    Tools::average( *field_manager.field(), field_manager.comm(), averages );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = (d+1) * random_numbers[1] * random_numbers[2]
                    * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( averages[d], local_val ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FieldTools, tensor_field_test )
{
    using namespace DataTransferKit;

    // Setup random numbers for the test.
    int num_rand = 4;
    double rand_max = 10.0;
    Teuchos::Array<double> random_numbers( num_rand );
    std::srand( 1 );
    for ( int i = 0; i < num_rand; ++i )
    {
        random_numbers[i] = rand_max * (double) std::rand() / RAND_MAX;
    }

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup a field manager.
    int field_dim = 9;
    typedef FieldTraits<ArrayField> FT;
    Teuchos::RCP<ArrayField> array_field =
        Teuchos::rcp( new ArrayField( field_dim*(my_rank+1), field_dim ) );
    FieldManager<ArrayField> field_manager( array_field, comm );

    // Test the field tools.
    typedef FieldTools<ArrayField> Tools;

    // Dimension iterators.
    FT::iterator dim_begin, dim_end;
    for ( int d = 0; d < field_dim; ++d )
    {
        dim_begin = Tools::dimBegin( *field_manager.field(), d );
        dim_end = Tools::dimEnd( *field_manager.field(), d );
        TEST_ASSERT( std::distance( dim_begin, dim_end ) == my_rank+1 );
    }

    FT::const_iterator const_dim_begin, const_dim_end;
    for ( int d = 0; d < field_dim; ++d )
    {
        const_dim_begin = Tools::dimBegin( *field_manager.field(), d );
        const_dim_end = Tools::dimEnd( *field_manager.field(), d );
        TEST_ASSERT( std::distance( const_dim_begin, const_dim_end )
                     == my_rank+1 );
    }

    // Length.
    int global_size = 0;
    for ( int i = 0; i < my_size; ++i )
    {
        global_size += field_dim*(i+1);
    }
    TEST_ASSERT(
        Tools::globalSize( *field_manager.field(), field_manager.comm() )
        == global_size );

    // Scalar filling.
    Tools::putScalar( *field_manager.field(), random_numbers[0] );

    // Views.
    double local_val;
    Teuchos::ArrayRCP<const double> const_view =
        Tools::view( *field_manager.field() );
    Teuchos::ArrayRCP<const double>::const_iterator const_view_iterator;
    for ( const_view_iterator = const_view.begin();
          const_view_iterator != const_view.end();
          ++const_view_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *const_view_iterator == local_val );
    }

    Teuchos::ArrayRCP<double> view =
        Tools::nonConstView( *field_manager.field() );
    Teuchos::ArrayRCP<double>::iterator view_iterator;
    for ( view_iterator = view.begin();
          view_iterator != view.end();
          ++view_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *view_iterator == local_val );
    }

    Teuchos::ArrayRCP<double> copy =
        Tools::copy( *field_manager.field() );
    Teuchos::ArrayRCP<double>::const_iterator copy_iterator;
    for ( copy_iterator = copy.begin();
          copy_iterator != copy.end();
          ++copy_iterator )
    {
        local_val = random_numbers[0];
        TEST_ASSERT( *copy_iterator == local_val );
    }

    // Multiple scalar filling.
    Teuchos::Array<double> fillers( field_dim );
    for ( int d = 0; d < field_dim; ++d )
    {
        fillers[d] = (d+1)*random_numbers[1];
    }
    Tools::putScalar( *field_manager.field(), fillers() );

    // Dim iterators.
    FT::const_iterator const_dim_iterator;
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( const_dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              const_dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++const_dim_iterator )
        {
            local_val = (d+1) * random_numbers[1];
            TEST_ASSERT( *const_dim_iterator == local_val );
        }
    }

    FT::iterator dim_iterator;
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1];
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Scaling.
    Tools::scale( *field_manager.field(), random_numbers[2] );
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1] * random_numbers[2];
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Multiple scaling.
    Teuchos::Array<double> multipliers( field_dim );
    for ( int d = 0; d < field_dim; ++d )
    {
        multipliers[d] = -(d+1) - random_numbers[3];
    }
    Tools::scale( *field_manager.field(), multipliers() );
    for ( int d = 0; d < field_dim; ++d )
    {
        for ( dim_iterator = Tools::dimBegin( *field_manager.field(), d );
              dim_iterator != Tools::dimEnd( *field_manager.field(), d );
              ++dim_iterator )
        {
            local_val = (d+1) * random_numbers[1]
                        * random_numbers[2] * (-(d+1) - random_numbers[3]);
            TEST_ASSERT( *dim_iterator == local_val );
        }
    }

    // Norms.
    Teuchos::Array<double> inf_norms;
    Tools::normInf( *field_manager.field(), field_manager.comm(), inf_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = (d+1) * random_numbers[1] * random_numbers[2]
                    * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( inf_norms[d], std::abs( local_val ) ) );
    }

    Teuchos::Array<double> one_norms;
    Tools::norm1( *field_manager.field(), field_manager.comm(), one_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                   (d+1) * random_numbers[1]
                   * random_numbers[2] * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( one_norms[d], std::abs( local_val ) ) );
    }

    Teuchos::Array<double> two_norms;
    Tools::norm2( *field_manager.field(), field_manager.comm(), two_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                   std::pow(
                       std::abs( (d+1) * random_numbers[1] * random_numbers[2]
                                 * (-(d+1) - random_numbers[3]) ), 2.0 );
        local_val = std::pow( local_val, 1.0/2.0 );
        TEST_ASSERT( softEquivalence( two_norms[d], local_val ) );
    }

    Teuchos::Array<double> q_norms;
    int q = std::floor( (10 * std::rand()) / RAND_MAX ) + 1;
    Tools::normQ( *field_manager.field(), field_manager.comm(), q, q_norms );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = ( global_size / FT::dim( *field_manager.field() ) ) *
                    std::pow(
                        std::abs( (d+1) * random_numbers[1] * random_numbers[2]
                                  * (-(d+1) - random_numbers[3]) ), q );
        local_val = std::pow( local_val, 1.0/q );
        TEST_ASSERT( softEquivalence( q_norms[d], local_val ) );
    }

    // Average.
    Teuchos::Array<double> averages;
    Tools::average( *field_manager.field(), field_manager.comm(), averages );
    for ( int d = 0; d < field_dim; ++d )
    {
        local_val = (d+1) * random_numbers[1] * random_numbers[2]
                    * (-(d+1) - random_numbers[3]);
        TEST_ASSERT( softEquivalence( averages[d], local_val ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FieldTools, 1d_coordinate_field_test )
{
    using namespace DataTransferKit;

    double double_limit = Teuchos::ScalarTraits<double>::rmax();

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();

    // Random number setup.
    int num_rand = 1000;
    double rand_max = 10.0;
    std::srand( my_rank*num_rand );

    // Setup a field manager.
    int field_dim = 1;
    typedef FieldTraits<ArrayField> FT;
    Teuchos::RCP<ArrayField> array_field =
        Teuchos::rcp( new ArrayField( (my_rank+1)*num_rand*field_dim, field_dim ) );
    FieldManager<ArrayField> field_manager( array_field, comm );

    // Test the field tools.
    typedef FieldTools<ArrayField> Tools;

    // Fill the field with random coordinates.
    FT::iterator field_iterator;
    for ( field_iterator = FT::begin( *field_manager.field() );
          field_iterator != FT::end( *field_manager.field() );
          ++field_iterator )
    {
        *field_iterator = rand_max * (double) std::rand() / RAND_MAX;
    }

    // Local bounding box.
    BoundingBox local_box =
        Tools::coordLocalBoundingBox( *field_manager.field() );
    double local_x_min =
        *std::min_element( Tools::dimBegin( *field_manager.field(), 0 ),
                           Tools::dimEnd( *field_manager.field(), 0 ) );
    double local_y_min = -double_limit;
    double local_z_min = -double_limit;
    double local_x_max =
        *std::max_element( Tools::dimBegin( *field_manager.field(), 0 ),
                           Tools::dimEnd( *field_manager.field(), 0 ) );
    double local_y_max = double_limit;
    double local_z_max = double_limit;

    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == local_x_min );
    TEST_ASSERT( local_bounds[1] == local_y_min );
    TEST_ASSERT( local_bounds[2] == local_z_min );
    TEST_ASSERT( local_bounds[3] == local_x_max );
    TEST_ASSERT( local_bounds[4] == local_y_max );
    TEST_ASSERT( local_bounds[5] == local_z_max );

    // Global bounding box.
    double global_x_min = 0;
    double global_x_max = 0;
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MIN,
                                    local_x_min,
                                    Teuchos::Ptr<double>(&global_x_min) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MAX,
                                    local_x_max,
                                    Teuchos::Ptr<double>(&global_x_max) );

    double global_y_min = -double_limit;
    double global_z_min = -double_limit;
    double global_y_max = double_limit;
    double global_z_max = double_limit;

    BoundingBox global_box =
        Tools::coordGlobalBoundingBox( *field_manager.field(),
                                       field_manager.comm() );
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == global_x_min );
    TEST_ASSERT( global_bounds[1] == global_y_min );
    TEST_ASSERT( global_bounds[2] == global_z_min );
    TEST_ASSERT( global_bounds[3] == global_x_max );
    TEST_ASSERT( global_bounds[4] == global_y_max );
    TEST_ASSERT( global_bounds[5] == global_z_max );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FieldTools, 2d_coordinate_field_test )
{
    using namespace DataTransferKit;

    double double_limit = Teuchos::ScalarTraits<double>::rmax();

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();

    // Random number setup.
    int num_rand = 1000;
    double rand_max = 10.0;
    std::srand( my_rank*num_rand );

    // Setup a field manager.
    int field_dim = 2;
    typedef FieldTraits<ArrayField> FT;
    Teuchos::RCP<ArrayField> array_field =
        Teuchos::rcp( new ArrayField( (my_rank+1)*num_rand*field_dim, field_dim ) );
    FieldManager<ArrayField> field_manager( array_field, comm );

    // Test the field tools.
    typedef FieldTools<ArrayField> Tools;

    // Fill the field with random coordinates.
    FT::iterator field_iterator;
    for ( field_iterator = FT::begin( *field_manager.field() );
          field_iterator != FT::end( *field_manager.field() );
          ++field_iterator )
    {
        *field_iterator = rand_max * (double) std::rand() / RAND_MAX;
    }

    // Local bounding box.
    BoundingBox local_box =
        Tools::coordLocalBoundingBox( *field_manager.field() );
    double local_x_min =
        *std::min_element( Tools::dimBegin( *field_manager.field(), 0 ),
                           Tools::dimEnd( *field_manager.field(), 0 ) );
    double local_y_min =
        *std::min_element( Tools::dimBegin( *field_manager.field(), 1 ),
                           Tools::dimEnd( *field_manager.field(), 1 ) );
    double local_z_min = -double_limit;
    double local_x_max =
        *std::max_element( Tools::dimBegin( *field_manager.field(), 0 ),
                           Tools::dimEnd( *field_manager.field(), 0 ) );
    double local_y_max =
        *std::max_element( Tools::dimBegin( *field_manager.field(), 1 ),
                           Tools::dimEnd( *field_manager.field(), 1 ) );
    double local_z_max = double_limit;

    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == local_x_min );
    TEST_ASSERT( local_bounds[1] == local_y_min );
    TEST_ASSERT( local_bounds[2] == local_z_min );
    TEST_ASSERT( local_bounds[3] == local_x_max );
    TEST_ASSERT( local_bounds[4] == local_y_max );
    TEST_ASSERT( local_bounds[5] == local_z_max );

    // Global bounding box.
    double global_x_min = 0;
    double global_y_min = 0;
    double global_x_max = 0;
    double global_y_max = 0;
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MIN,
                                    local_x_min,
                                    Teuchos::Ptr<double>(&global_x_min) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MIN,
                                    local_y_min,
                                    Teuchos::Ptr<double>(&global_y_min) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MAX,
                                    local_x_max,
                                    Teuchos::Ptr<double>(&global_x_max) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MAX,
                                    local_y_max,
                                    Teuchos::Ptr<double>(&global_y_max) );

    double global_z_min = -double_limit;
    double global_z_max = double_limit;

    BoundingBox global_box =
        Tools::coordGlobalBoundingBox( *field_manager.field(),
                                       field_manager.comm() );
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == global_x_min );
    TEST_ASSERT( global_bounds[1] == global_y_min );
    TEST_ASSERT( global_bounds[2] == global_z_min );
    TEST_ASSERT( global_bounds[3] == global_x_max );
    TEST_ASSERT( global_bounds[4] == global_y_max );
    TEST_ASSERT( global_bounds[5] == global_z_max );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( FieldTools, 3d_coordinate_field_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();

    // Random number setup.
    int num_rand = 1000;
    double rand_max = 10.0;
    std::srand( my_rank*num_rand );

    // Setup a field manager.
    int field_dim = 3;
    typedef FieldTraits<ArrayField> FT;
    Teuchos::RCP<ArrayField> array_field =
        Teuchos::rcp( new ArrayField( (my_rank+1)*num_rand*field_dim, field_dim ) );
    FieldManager<ArrayField> field_manager( array_field, comm );

    // Test the field tools.
    typedef FieldTools<ArrayField> Tools;

    // Fill the field with random coordinates.
    FT::iterator field_iterator;
    for ( field_iterator = FT::begin( *field_manager.field() );
          field_iterator != FT::end( *field_manager.field() );
          ++field_iterator )
    {
        *field_iterator = rand_max * (double) std::rand() / RAND_MAX;
    }

    // Local bounding box.
    BoundingBox local_box =
        Tools::coordLocalBoundingBox( *field_manager.field() );
    double local_x_min =
        *std::min_element( Tools::dimBegin( *field_manager.field(), 0 ),
                           Tools::dimEnd( *field_manager.field(), 0 ) );
    double local_y_min =
        *std::min_element( Tools::dimBegin( *field_manager.field(), 1 ),
                           Tools::dimEnd( *field_manager.field(), 1 ) );
    double local_z_min =
        *std::min_element( Tools::dimBegin( *field_manager.field(), 2 ),
                           Tools::dimEnd( *field_manager.field(), 2 ) );
    double local_x_max =
        *std::max_element( Tools::dimBegin( *field_manager.field(), 0 ),
                           Tools::dimEnd( *field_manager.field(), 0 ) );
    double local_y_max =
        *std::max_element( Tools::dimBegin( *field_manager.field(), 1 ),
                           Tools::dimEnd( *field_manager.field(), 1 ) );
    double local_z_max =
        *std::max_element( Tools::dimBegin( *field_manager.field(), 2 ),
                           Tools::dimEnd( *field_manager.field(), 2 ) );

    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == local_x_min );
    TEST_ASSERT( local_bounds[1] == local_y_min );
    TEST_ASSERT( local_bounds[2] == local_z_min );
    TEST_ASSERT( local_bounds[3] == local_x_max );
    TEST_ASSERT( local_bounds[4] == local_y_max );
    TEST_ASSERT( local_bounds[5] == local_z_max );

    // Global bounding box.
    double global_x_min = 0;
    double global_y_min = 0;
    double global_z_min = 0;
    double global_x_max = 0;
    double global_y_max = 0;
    double global_z_max = 0;
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MIN,
                                    local_x_min,
                                    Teuchos::Ptr<double>(&global_x_min) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MIN,
                                    local_y_min,
                                    Teuchos::Ptr<double>(&global_y_min) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MIN,
                                    local_z_min,
                                    Teuchos::Ptr<double>(&global_z_min) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MAX,
                                    local_x_max,
                                    Teuchos::Ptr<double>(&global_x_max) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MAX,
                                    local_y_max,
                                    Teuchos::Ptr<double>(&global_y_max) );
    Teuchos::reduceAll<int,double>( *comm, Teuchos::REDUCE_MAX,
                                    local_z_max,
                                    Teuchos::Ptr<double>(&global_z_max) );

    BoundingBox global_box =
        Tools::coordGlobalBoundingBox( *field_manager.field(),
                                       field_manager.comm() );
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == global_x_min );
    TEST_ASSERT( global_bounds[1] == global_y_min );
    TEST_ASSERT( global_bounds[2] == global_z_min );
    TEST_ASSERT( global_bounds[3] == global_x_max );
    TEST_ASSERT( global_bounds[4] == global_y_max );
    TEST_ASSERT( global_bounds[5] == global_z_max );
}

//---------------------------------------------------------------------------//
// end tstFieldTools.cpp
//---------------------------------------------------------------------------//

