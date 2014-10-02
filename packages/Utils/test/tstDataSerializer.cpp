//---------------------------------------------------------------------------//
/*!
 * \file   tstSerializer.cpp
 * \author Stuart Slattery
 * \brief  Serializer class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_DataSerializer.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_as.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
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
// Basic Struct.
//---------------------------------------------------------------------------//
struct DataHolder
{
    double data;
};

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( DataSerializer, size_test )
{
    bool data_bool = true;
    unsigned int data_uint = 1;
    int data_int = -4;
    float data_flt = -0.4332;
    double data_dbl = 3.2;
    DataHolder data_holder;

    std::size_t buffer_size = sizeof(bool) + sizeof(unsigned int) +
			      sizeof(int) + sizeof(float) +
			      sizeof(double) + sizeof(DataHolder);

    DataTransferKit::DataSerializer serializer;
    serializer.computeBufferSizeMode();
    serializer.pack( data_bool );
    serializer.pack( data_uint );
    serializer.pack( data_int );
    serializer.pack( data_flt );
    serializer.pack( data_dbl );
    serializer.pack( data_holder );

    TEST_EQUALITY( serializer.size(), buffer_size );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( DataSerializer, stream_size_test )
{
    bool data_bool = true;
    unsigned int data_uint = 1;
    int data_int = -4;
    float data_flt = -0.4332;
    double data_dbl = 3.2;
    DataHolder data_holder;

    std::size_t buffer_size = sizeof(bool) + sizeof(unsigned int) +
			      sizeof(int) + sizeof(float) +
			      sizeof(double) + sizeof(DataHolder);

    DataTransferKit::DataSerializer serializer;
    serializer.computeBufferSizeMode();
    serializer << data_bool << data_uint << data_int << data_flt 
	       << data_dbl << data_holder;

    TEST_EQUALITY( serializer.size(), buffer_size );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( DataSerializer, pack_unpack_test )
{
    bool data_bool = true;
    unsigned int data_uint = 1;
    int data_int = -4;
    float data_flt = -0.4332;
    double data_dbl = 3.2;
    DataHolder data_holder;
    double data_val = 2;
    data_holder.data = data_val;

    DataTransferKit::DataSerializer serializer;
    serializer.computeBufferSizeMode();
    serializer << data_bool << data_uint << data_int << data_flt 
	       << data_dbl << data_holder;

    Teuchos::Array<char> buffer( serializer.size() );
    serializer.setBuffer( buffer.size(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.getPtr(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.begin(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.end(), buffer.getRawPtr()+serializer.size() );

    serializer << data_bool << data_uint << data_int << data_flt 
	       << data_dbl << data_holder;

    DataTransferKit::DataDeserializer deserializer;
    deserializer.setBuffer( buffer.size(), buffer.getRawPtr() );

    bool ds_bool = false;
    unsigned int ds_uint = 0;
    int ds_int = 0;
    float ds_flt = 0.0;
    double ds_dbl = 0.0;
    DataHolder ds_holder;

    deserializer.unpack( ds_bool );
    TEST_EQUALITY( ds_bool, data_bool );

    deserializer.unpack( ds_uint );
    TEST_EQUALITY( ds_uint, data_uint );

    deserializer.unpack( ds_int );
    TEST_EQUALITY( ds_int, data_int );

    deserializer.unpack( ds_flt );
    TEST_EQUALITY( ds_flt, data_flt );

    deserializer.unpack( ds_dbl );
    TEST_EQUALITY( ds_dbl, data_dbl );

    deserializer.unpack( ds_holder );
    TEST_EQUALITY( ds_holder.data, data_val );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( DataSerializer, stream_pack_unpack_test )
{
    bool data_bool = true;
    unsigned int data_uint = 1;
    int data_int = -4;
    float data_flt = -0.4332;
    double data_dbl = 3.2;
    DataHolder data_holder;
    double data_val = 2;
    data_holder.data = data_val;

    DataTransferKit::DataSerializer serializer;
    serializer.computeBufferSizeMode();
    serializer << data_bool << data_uint << data_int << data_flt 
	       << data_dbl << data_holder;

    Teuchos::Array<char> buffer( serializer.size() );
    serializer.setBuffer( buffer.size(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.getPtr(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.begin(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.end(), buffer.getRawPtr()+serializer.size() );

    serializer << data_bool << data_uint << data_int << data_flt 
	       << data_dbl << data_holder;

    DataTransferKit::DataDeserializer deserializer;
    deserializer.setBuffer( buffer.size(), buffer.getRawPtr() );

    bool ds_bool = false;
    unsigned int ds_uint = 0;
    int ds_int = 0;
    float ds_flt = 0.0;
    double ds_dbl = 0.0;
    DataHolder ds_holder;

    deserializer >> ds_bool;
    TEST_EQUALITY( ds_bool, data_bool );

    deserializer >> ds_uint;
    TEST_EQUALITY( ds_uint, data_uint );

    deserializer >> ds_int;
    TEST_EQUALITY( ds_int, data_int );

    deserializer >> ds_flt;
    TEST_EQUALITY( ds_flt, data_flt );

    deserializer >> ds_dbl;
    TEST_EQUALITY( ds_dbl, data_dbl );

    deserializer >> ds_holder;
    TEST_EQUALITY( ds_holder.data, data_val );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( DataSerializer, view_pack_unpack_test )
{
    bool data_bool = true;
    unsigned int data_uint = 1;
    int data_int = -4;
    float data_flt = -0.4332;
    double data_dbl = 3.2;
    DataHolder data_holder;
    double data_val = 2;
    data_holder.data = data_val;

    DataTransferKit::DataSerializer serializer;
    serializer.computeBufferSizeMode();
    serializer << data_bool << data_uint << data_int << data_flt 
	       << data_dbl << data_holder;

    Teuchos::Array<char> buffer( serializer.size() );
    serializer.setBuffer( buffer() );
    TEST_EQUALITY( serializer.getPtr(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.begin(), buffer.getRawPtr() );
    TEST_EQUALITY( serializer.end(), buffer.getRawPtr()+serializer.size() );

    serializer << data_bool << data_uint << data_int << data_flt 
	       << data_dbl << data_holder;

    DataTransferKit::DataDeserializer deserializer;
    deserializer.setBuffer( buffer() );

    bool ds_bool = false;
    unsigned int ds_uint = 0;
    int ds_int = 0;
    float ds_flt = 0.0;
    double ds_dbl = 0.0;
    DataHolder ds_holder;

    deserializer >> ds_bool;
    TEST_EQUALITY( ds_bool, data_bool );

    deserializer >> ds_uint;
    TEST_EQUALITY( ds_uint, data_uint );

    deserializer >> ds_int;
    TEST_EQUALITY( ds_int, data_int );

    deserializer >> ds_flt;
    TEST_EQUALITY( ds_flt, data_flt );

    deserializer >> ds_dbl;
    TEST_EQUALITY( ds_dbl, data_dbl );

    deserializer >> ds_holder;
    TEST_EQUALITY( ds_holder.data, data_val );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( DataSerializer, broadcast_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();

    bool data_bool = true;
    unsigned int data_uint = 1;
    int data_int = -4;
    float data_flt = -0.4332;
    double data_dbl = 3.2;
    DataHolder data_holder;
    double data_val = 2;
    data_holder.data = data_val;

    std::size_t buffer_size = sizeof(bool) + sizeof(unsigned int) +
			      sizeof(int) + sizeof(float) +
			      sizeof(double) + sizeof(DataHolder);

    Teuchos::Array<char> buffer( buffer_size );

    if ( comm_rank == 0 )
    {
	DataTransferKit::DataSerializer serializer;
	serializer.setBuffer( buffer.size(), buffer.getRawPtr() );
	serializer << data_bool << data_uint << data_int << data_flt 
		   << data_dbl << data_holder;
    }
    comm->barrier();

    Teuchos::broadcast( *comm, 0, buffer() );

    DataTransferKit::DataDeserializer deserializer;
    deserializer.setBuffer( buffer.size(), buffer.getRawPtr() );

    bool ds_bool = false;
    unsigned int ds_uint = 0;
    int ds_int = 0;
    float ds_flt = 0.0;
    double ds_dbl = 0.0;
    DataHolder ds_holder;

    deserializer >> ds_bool;
    TEST_EQUALITY( ds_bool, data_bool );

    deserializer >> ds_uint;
    TEST_EQUALITY( ds_uint, data_uint );

    deserializer >> ds_int;
    TEST_EQUALITY( ds_int, data_int );

    deserializer >> ds_flt;
    TEST_EQUALITY( ds_flt, data_flt );

    deserializer >> ds_dbl;
    TEST_EQUALITY( ds_dbl, data_dbl );

    deserializer >> ds_holder;
    TEST_EQUALITY( ds_holder.data, data_val );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( DataSerializer, view_broadcast_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();

    bool data_bool = true;
    unsigned int data_uint = 1;
    int data_int = -4;
    float data_flt = -0.4332;
    double data_dbl = 3.2;
    DataHolder data_holder;
    double data_val = 2;
    data_holder.data = data_val;

    std::size_t buffer_size = sizeof(bool) + sizeof(unsigned int) +
			      sizeof(int) + sizeof(float) +
			      sizeof(double) + sizeof(DataHolder);

    Teuchos::Array<char> buffer( buffer_size );

    if ( comm_rank == 0 )
    {
	DataTransferKit::DataSerializer serializer;
	serializer.setBuffer( buffer() );
	serializer << data_bool << data_uint << data_int << data_flt 
		   << data_dbl << data_holder;
    }
    comm->barrier();

    Teuchos::broadcast( *comm, 0, buffer() );

    DataTransferKit::DataDeserializer deserializer;
    deserializer.setBuffer( buffer() );

    bool ds_bool = false;
    unsigned int ds_uint = 0;
    int ds_int = 0;
    float ds_flt = 0.0;
    double ds_dbl = 0.0;
    DataHolder ds_holder;

    deserializer >> ds_bool;
    TEST_EQUALITY( ds_bool, data_bool );

    deserializer >> ds_uint;
    TEST_EQUALITY( ds_uint, data_uint );

    deserializer >> ds_int;
    TEST_EQUALITY( ds_int, data_int );

    deserializer >> ds_flt;
    TEST_EQUALITY( ds_flt, data_flt );

    deserializer >> ds_dbl;
    TEST_EQUALITY( ds_dbl, data_dbl );

    deserializer >> ds_holder;
    TEST_EQUALITY( ds_holder.data, data_val );
}

//---------------------------------------------------------------------------//
// end tstDataSerializer.cpp
//---------------------------------------------------------------------------//
