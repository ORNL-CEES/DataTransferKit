//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstDataField.cpp
 * \author Stuart Slattery
 * \date   Fri Nov 18 14:43:10 2011
 * \brief  DataField unit tests
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Mesh_Point.hpp>

#include <Coupler_DataSource.hpp>
#include <Coupler_DataTarget.hpp>
#include <Coupler_DataField.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
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
// DATA CLASS
//---------------------------------------------------------------------------//
class Data_Container
{
  public:

    typedef Coupler::Point<int,double>     PointType;

  private:

    double scalar_data;
    std::vector<double> distributed_data;
    std::vector<PointType> distributed_points;

  public:

    Data_Container()
    { /* ... */ }

    ~Data_Container()
    { /* ... */ }

    void set_distributed_points(std::vector<PointType> points)
    {
	distributed_points = points;
    }

    Teuchos::ArrayView<PointType> get_distributed_points()
    {
	return Teuchos::ArrayView<PointType>(distributed_points);
    }

    Teuchos::ArrayView<double> get_distributed_data()
    {
	distributed_data.resize(5);
	return Teuchos::ArrayView<double>(distributed_data);
    }

    void set_scalar_data(double data)
    {
	scalar_data = data;
    }

    double get_scalar_data()
    {
	return scalar_data;
    }
};

//---------------------------------------------------------------------------//
// INTERFACE IMPLEMENTATIONS
//---------------------------------------------------------------------------//

namespace Coupler {

// transfer data source implementation - this implementation specifies double
// as the data type
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class test_DataSource 
    : public DataSource<DataType_T, HandleType_T, CoordinateType_T>
{
  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef Point<HandleType,CoordinateType>         PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;

  private:

    std::vector<double> private_data;
    Teuchos::RCP<Data_Container> container;
    std::vector<PointType> local_points;
    int myRank;
    int mySize;

  public:

    test_DataSource(Teuchos::RCP<Data_Container> _container)
	: container(_container)
    { 
	myRank = getDefaultComm<OrdinalType>()->getRank();
	mySize = getDefaultComm<OrdinalType>()->getSize();
    }

    ~test_DataSource()
    { /* ... */ }

    RCP_Communicator comm()
    {
	return getDefaultComm<OrdinalType>();
    }

    bool field_supported(const std::string &field_name)
    {
	bool return_val = false;

	if (field_name == "DISTRIBUTED_TEST_FIELD")
	{
	    return_val = true;
	}

	else if (field_name == "SCALAR_TEST_FIELD")
	{
	    return_val = true;
	}

	return return_val;
    }

    bool get_points(const PointType &point)
    {
	bool return_val = false;

	if ( point.x() == 1.0*(mySize-myRank-1) && 
	     point.y() == 2.0*(mySize-myRank-1) && 
	     point.z() == 3.0*(mySize-myRank-1) )
	{
	    return_val = true;
	    local_points.push_back(point);
	    container->set_distributed_points(local_points);
	}

	return return_val;
    }

    const Teuchos::ArrayView<double> send_data(const std::string &field_name)
    {
	Teuchos::ArrayView<double> return_view;

	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    private_data.resize(local_points.size());
	    std::fill( private_data.begin(), private_data.end(), 1.0*myRank);
	    Teuchos::ArrayView<double> private_view(private_data);
	    return_view =  private_view;
	}

	return return_view;
    }

    double set_global_data(const std::string &field_name)
    {
	double return_val = 0.0;

	if ( field_name == "SCALAR_TEST_FIELD" )
	{
	    return_val = 5.352;
	}

	return return_val;
    }
};

//---------------------------------------------------------------------------//
// transfer data target implementation - this implementation specifies double
// as the data type
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class test_DataTarget 
    : public DataTarget<DataType_T, HandleType_T, CoordinateType_T>
{
  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef Point<HandleType,CoordinateType>         PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;

  private:

    Teuchos::RCP<Data_Container> container;
    std::vector<PointType> local_points;
    int myRank;
    int mySize;

  public:

    test_DataTarget(Teuchos::RCP<Data_Container> _container)
	: container(_container)
    { 
	myRank = getDefaultComm<OrdinalType>()->getRank();
	mySize = getDefaultComm<OrdinalType>()->getSize();
    }

    ~test_DataTarget()
    { /* ... */ }

    RCP_Communicator comm()
    {
	return getDefaultComm<OrdinalType>();
    }

    bool field_supported(const std::string &field_name)
    {
	bool return_val = false;

	if (field_name == "DISTRIBUTED_TEST_FIELD")
	{
	    return_val = true;
	}

	else if (field_name == "SCALAR_TEST_FIELD")
	{
	    return_val = true;
	}

	return return_val;
    }

    const Teuchos::ArrayView<PointType> 
    set_points(const std::string &field_name)
    {
	Teuchos::ArrayView<PointType> return_view;

	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    for (int i = 0; i < 5; ++i)
	    {
		int this_handle = 5*myRank + i;
		PointType point(this_handle, 1.0*myRank, 
				2.0*myRank, 3.0*myRank);
		local_points.push_back(point);
	    }
	    return_view = Teuchos::ArrayView<PointType>(local_points);
	}

	return return_view;
    }

    Teuchos::ArrayView<DataType> receive_data(const std::string &field_name)
    {
	Teuchos::ArrayView<DataType> return_view;

	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    return_view = container->get_distributed_data();
	}

	return return_view;
    }

    void get_global_data(const std::string &field_name,
			 const DataType &data)
    {
	if ( field_name == "SCALAR_TEST_FIELD" )
	{
	    container->set_scalar_data(data);
	}
    }
};

} // end namespace Coupler

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace Coupler {

TEUCHOS_UNIT_TEST( DataField, distributed_container_test )
{
    // create a data container instance for checking the data under the source
    // interface.
    Teuchos::RCP<Data_Container> source_container 
	= Teuchos::rcp(new Data_Container);

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> target_container 
	= Teuchos::rcp(new Data_Container);

    // Create an instance of the source interface.
    Teuchos::RCP<DataSource<double,int,double> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double>(source_container));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double>(target_container));

    // Create a distributed field for these interfaces to be transferred.
    DataField<double,int,double> field(getDefaultComm<int>(),
				       "DISTRIBUTED_TEST_FIELD", 
				       "DISTRIBUTED_TEST_FIELD",
				       tds, 
				       tdt);

    // Test the basic container functionality of the field.
    TEST_ASSERT( field.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( field.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( field.source_name() == "DISTRIBUTED_TEST_FIELD" );
    TEST_ASSERT( field.source() == tds );
    TEST_ASSERT( field.target_name() == "DISTRIBUTED_TEST_FIELD" );
    TEST_ASSERT( field.target() == tdt );
    TEST_ASSERT( !field.is_scalar() );
    TEST_ASSERT( field.is_mapped() );
}

TEUCHOS_UNIT_TEST( DataField, scalar_container_test )
{
    // create a data container instance for checking the data under the source
    // interface.
    Teuchos::RCP<Data_Container> source_container 
	= Teuchos::rcp(new Data_Container);

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> target_container 
	= Teuchos::rcp(new Data_Container);

    // Create an instance of the source interface.
    Teuchos::RCP<DataSource<double,int,double> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double>(source_container));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double>(target_container));

    // Create a distributed field for these interfaces to be transferred.
    DataField<double,int,double> field(getDefaultComm<int>(),
				       "SCALAR_TEST_FIELD", 
				       "SCALAR_TEST_FIELD", 
				       tds, 
				       tdt,
				       true);

    // Test the basic container functionality of the field.
    TEST_ASSERT( field.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( field.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( field.source_name() == "SCALAR_TEST_FIELD" );
    TEST_ASSERT( field.source() == tds );
    TEST_ASSERT( field.target_name() == "SCALAR_TEST_FIELD" );
    TEST_ASSERT( field.target() == tdt );
    TEST_ASSERT( field.is_scalar() );
}

TEUCHOS_UNIT_TEST( DataField, mapping_test )
{
    // create a data container instance for checking the data under the source
    // interface.
    Teuchos::RCP<Data_Container> source_container 
	= Teuchos::rcp(new Data_Container);

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> target_container 
	= Teuchos::rcp(new Data_Container);

    // Create an instance of the source interface.
    Teuchos::RCP<DataSource<double,int,double> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double>(source_container));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double>(target_container));

    // Create a distributed field for these interfaces to be transferred.
    DataField<double,int,double> field(getDefaultComm<int>(),
				       "DISTRIBUTED_TEST_FIELD", 
				       "DISTRIBUTED_TEST_FIELD",
				       tds, 
				       tdt);

    // Check the points under the source interface to the verify communication
    // pattern of sending target points to the source.
    int myRank = getDefaultComm<int>()->getRank();
    int mySize = getDefaultComm<int>()->getSize();
    int flippedRank = mySize-myRank-1;
    TEST_ASSERT( field.is_mapped() );
    TEST_ASSERT( source_container->get_distributed_points().size() == 5 );
    for (int i = 0; i < 5; ++i)
    {
	TEST_ASSERT( source_container->get_distributed_points()[i].handle() 
		     == 5*flippedRank+i );
	TEST_ASSERT( source_container->get_distributed_points()[i].x()
		     == 1.0*flippedRank );
	TEST_ASSERT( source_container->get_distributed_points()[i].y()
		     == 2.0*flippedRank );
	TEST_ASSERT( source_container->get_distributed_points()[i].z()
		     == 3.0*flippedRank );
    }

    // Check the Tpetra maps for both the source and the target.
    TEST_ASSERT( (int) field.source_map()->getGlobalNumElements() 
		 == mySize*5 );
    TEST_ASSERT( (int) field.source_map()->getNodeNumElements() == 5 );
    for (int i = 0; i < 5; ++i)
    {
	TEST_ASSERT( field.source_map()->getNodeElementList()[i] 
		     == 5*flippedRank+i );
    }

    TEST_ASSERT( (int) field.target_map()->getGlobalNumElements() 
		 == mySize*5 );
    TEST_ASSERT( (int) field.target_map()->getNodeNumElements() == 5 );
    for (int i = 0; i < 5; ++i)
    {
	TEST_ASSERT( field.target_map()->getNodeElementList()[i] 
		     == 5*myRank+i );
    }
}

TEUCHOS_UNIT_TEST( DataField, Scalar_Transfer_Test )
{
    // create a data container instance for checking the data under the source
    // interface.
    Teuchos::RCP<Data_Container> source_container =
	Teuchos::rcp(new Data_Container);

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> target_container 
	= Teuchos::rcp(new Data_Container);

    // Create an instance of the source interface.
    Teuchos::RCP<DataSource<double,int,double> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double>(source_container));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double>(target_container));

    // Create a distributed field for these interfaces to be transferred.
    DataField<double,int,double> field(getDefaultComm<int>(),
				       "SCALAR_TEST_FIELD", 
				       "SCALAR_TEST_FIELD", 
				       tds, 
				       tdt,
				       true);

    // Do scalar transfer.
    field.transfer();
    TEST_ASSERT( target_container->get_scalar_data() == 5.352 );
}

TEUCHOS_UNIT_TEST( DataField, Distributed_Transfer_Test )
{
    // create a data container instance for checking the data under the source
    // interface.
    Teuchos::RCP<Data_Container> source_container 
	= Teuchos::rcp(new Data_Container);

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> target_container 
	= Teuchos::rcp(new Data_Container);

    // Create an instance of the source interface.
    Teuchos::RCP<DataSource<double,int,double> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double>(source_container));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double>(target_container));

    // Create a distributed field for these interfaces to be transferred.
    DataField<double,int,double> field(getDefaultComm<int>(),
				       "DISTRIBUTED_TEST_FIELD", 
				       "DISTRIBUTED_TEST_FIELD",
				       tds, 
				       tdt);

    // Do the transfer.
    field.transfer();

    // Check the transferred data under the target interface.
    TEST_ASSERT( target_container->get_distributed_data().size() == 5 );
    int myRank = getDefaultComm<int>()->getRank();
    int mySize = getDefaultComm<int>()->getSize();
    int flippedRank = mySize-myRank-1;
    for (int i = 0; i < 5; ++i)
    {
	TEST_ASSERT( target_container->get_distributed_data()[i]
		     == 1.0*flippedRank );
    }
}

//---------------------------------------------------------------------------//

} // end namespace Coupler

//---------------------------------------------------------------------------//
//                        end of tstDataField.cpp
//---------------------------------------------------------------------------//
