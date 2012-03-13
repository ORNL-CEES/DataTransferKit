//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstCopyOperator.cpp
 * \author Stuart Slattery
 * \date   Fri Nov 18 14:43:10 2011
 * \brief  CopyOperator unit tests
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Coupler_Point.hpp>
#include <Coupler_DataSource.hpp>
#include <Coupler_DataTarget.hpp>
#include <Coupler_CopyOperator.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ENull.hpp"
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

    typedef Coupler::Point<3,int,double>     PointType;

  private:

    double scalar_data;
    std::vector<double> distributed_data;
    std::vector<PointType> distributed_points;

  public:

    Data_Container()
	: scalar_data(0.0)
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
template<class DataType, class HandleType, class CoordinateType, int DIM>
class test_DataSource 
    : public DataSource<DataType,HandleType,CoordinateType,DIM>
{
  public:

    typedef int                                      OrdinalType;
    typedef Point<DIM,HandleType,CoordinateType>     PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;

  private:

    std::vector<double> private_data;
    Teuchos::RCP<Data_Container> container;
    std::vector<PointType> local_points;
    int myRank;
    int mySize;
    RCP_Communicator comm;

  public:

    test_DataSource(Teuchos::RCP<Data_Container> _container,
		    RCP_Communicator _comm)
	: container(_container)
	, comm(_comm)
    { 
	myRank = comm->getRank();
	mySize = comm->getSize();
    }

    ~test_DataSource()
    { /* ... */ }

    RCP_Communicator get_source_comm()
    {
	return comm;
    }

    bool is_field_supported(const std::string &field_name)
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

    bool is_local_point(const PointType &test_point)
    {
	bool return_val = false;

	Teuchos::Tuple<CoordinateType,DIM> test_coords = test_point.getCoords();
	if ( test_coords[0] == 1.0*(mySize-myRank-1) && 
	     test_coords[1] == 2.0*(mySize-myRank-1) && 
	     test_coords[2] == 3.0*(mySize-myRank-1) )
	{
	    return_val = true;
	    local_points.push_back(test_point);
	    container->set_distributed_points(local_points);
	}

	return return_val;
    }

    const Teuchos::ArrayView<double> 
    get_source_data(const std::string &field_name)
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

    double get_global_source_data(const std::string &field_name)
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
template<class DataType, class HandleType, class CoordinateType, int DIM>
class test_DataTarget 
    : public DataTarget<DataType,HandleType,CoordinateType,DIM>
{
  public:

    typedef int                                      OrdinalType;
    typedef Point<DIM,HandleType,CoordinateType>     PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;

  private:

    Teuchos::RCP<Data_Container> container;
    std::vector<PointType> local_points;
    int myRank;
    int mySize;
    RCP_Communicator comm;

  public:

    test_DataTarget( Teuchos::RCP<Data_Container> _container,
		     RCP_Communicator _comm )
	: container(_container)
	, comm(_comm)
    { 
	myRank = comm->getRank();
	mySize = comm->getSize();
    }

    ~test_DataTarget()
    { /* ... */ }

    RCP_Communicator get_target_comm()
    {
	return comm;
    }

    bool is_field_supported(const std::string &field_name)
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
    get_target_points(const std::string &field_name)
    {
	Teuchos::ArrayView<PointType> return_view;

	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    for (int i = 0; i < 5; ++i)
	    {
		int new_handle = 5*myRank + i;
		PointType new_point = point(new_handle, 
					    1.0*myRank, 
					    2.0*myRank, 
					    3.0*myRank);
		local_points.push_back(new_point);
	    }
	    return_view = Teuchos::ArrayView<PointType>(local_points);
	}

	return return_view;
    }

    Teuchos::ArrayView<DataType> 
    get_target_data_space(const std::string &field_name)
    {
	Teuchos::ArrayView<DataType> return_view;

	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    return_view = container->get_distributed_data();
	}

	return return_view;
    }

    void set_global_target_data(const std::string &field_name,
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

TEUCHOS_UNIT_TEST( CopyOperator, distributed_container_test )
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
    Teuchos::RCP<DataSource<double,int,double,3> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double,3>(source_container,
							      getDefaultComm<int>()));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double,3>(target_container,
							      getDefaultComm<int>()));

    // Create a distributed field for these interfaces to be transferred.
    CopyOperator<double,int,double,3> field_op( getDefaultComm<int>(),
						"DISTRIBUTED_TEST_FIELD", 
						"DISTRIBUTED_TEST_FIELD",
						tds, 
						tdt );

    // Test the basic container functionality of the field.
    TEST_ASSERT( field_op.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( field_op.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( field_op.source_name() == "DISTRIBUTED_TEST_FIELD" );
    TEST_ASSERT( field_op.source() == tds );
    TEST_ASSERT( field_op.target_name() == "DISTRIBUTED_TEST_FIELD" );
    TEST_ASSERT( field_op.target() == tdt );
    TEST_ASSERT( !field_op.is_scalar() );
    TEST_ASSERT( !field_op.is_mapped() );
}

TEUCHOS_UNIT_TEST( CopyOperator, scalar_container_test )
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
    Teuchos::RCP<DataSource<double,int,double,3> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double,3>(source_container,
							      getDefaultComm<int>()));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double,3>(target_container,
							      getDefaultComm<int>()));

    // Create a distributed field for these interfaces to be transferred.
    CopyOperator<double,int,double,3> field_op( getDefaultComm<int>(),
						"SCALAR_TEST_FIELD", 
						"SCALAR_TEST_FIELD", 
						tds, 
						tdt,
						true );

    // Test the basic container functionality of the field.
    TEST_ASSERT( field_op.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( field_op.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( field_op.source_name() == "SCALAR_TEST_FIELD" );
    TEST_ASSERT( field_op.source() == tds );
    TEST_ASSERT( field_op.target_name() == "SCALAR_TEST_FIELD" );
    TEST_ASSERT( field_op.target() == tdt );
    TEST_ASSERT( field_op.is_scalar() );
}

TEUCHOS_UNIT_TEST( CopyOperator, mapping_test )
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
    Teuchos::RCP<DataSource<double,int,double,3> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double,3>(source_container,
							      getDefaultComm<int>()));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double,3>(target_container,
							      getDefaultComm<int>()));

    // Create a distributed field for these interfaces to be transferred.
    CopyOperator<double,int,double,3> field_op( getDefaultComm<int>(),
						"DISTRIBUTED_TEST_FIELD", 
						"DISTRIBUTED_TEST_FIELD",
						tds, 
						tdt );

    // Generate the mapping.
    field_op.create_copy_mapping();

    // Check the points under the source interface to the verify communication
    // pattern of sending target points to the source.
    int myRank = getDefaultComm<int>()->getRank();
    int mySize = getDefaultComm<int>()->getSize();
    int flippedRank = mySize-myRank-1;
    TEST_ASSERT( field_op.is_mapped() );
    TEST_ASSERT( source_container->get_distributed_points().size() == 5 );
    for (int i = 0; i < 5; ++i)
    {
	TEST_ASSERT( source_container->get_distributed_points()[i].getHandle() 
		     == 5*flippedRank+i );
	TEST_ASSERT( source_container->get_distributed_points()[i].getCoords()[0]
		     == 1.0*flippedRank );
	TEST_ASSERT( source_container->get_distributed_points()[i].getCoords()[1]
		     == 2.0*flippedRank );
	TEST_ASSERT( source_container->get_distributed_points()[i].getCoords()[2]
		     == 3.0*flippedRank );
    }

    // Check the Tpetra maps for both the source and the target.
    TEST_ASSERT( (int) field_op.source_map()->getGlobalNumElements() 
		 == mySize*5 );
    TEST_ASSERT( (int) field_op.source_map()->getNodeNumElements() == 5 );
    for (int i = 0; i < 5; ++i)
    {
	TEST_ASSERT( field_op.source_map()->getNodeElementList()[i] 
		     == 5*flippedRank+i );
    }

    TEST_ASSERT( (int) field_op.target_map()->getGlobalNumElements() 
		 == mySize*5 );
    TEST_ASSERT( (int) field_op.target_map()->getNodeNumElements() == 5 );
    for (int i = 0; i < 5; ++i)
    {
	TEST_ASSERT( field_op.target_map()->getNodeElementList()[i] 
		     == 5*myRank+i );
    }
}

TEUCHOS_UNIT_TEST( CopyOperator, Scalar_Transfer_Test )
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
    Teuchos::RCP<DataSource<double,int,double,3> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double,3>(source_container,
							      getDefaultComm<int>()));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double,3>(target_container,
							      getDefaultComm<int>()));

    // Create a distributed field for these interfaces to be transferred.
    CopyOperator<double,int,double,3> field_op( getDefaultComm<int>(),
						"SCALAR_TEST_FIELD", 
						"SCALAR_TEST_FIELD", 
						tds, 
						tdt,
						true );

    // Do scalar transfer.
    field_op.copy();

    TEST_ASSERT( target_container->get_scalar_data() == 5.352 );
}

TEUCHOS_UNIT_TEST( CopyOperator, Distributed_Transfer_Test )
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
    Teuchos::RCP<DataSource<double,int,double,3> > tds = 
	Teuchos::rcp(new test_DataSource<double,int,double,3>(source_container,
							      getDefaultComm<int>()));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double,3>(target_container,
							      getDefaultComm<int>()));

    // Create a distributed field for these interfaces to be transferred.
    CopyOperator<double,int,double,3> field_op( getDefaultComm<int>(),
						"DISTRIBUTED_TEST_FIELD", 
						"DISTRIBUTED_TEST_FIELD",
						tds, 
						tdt);

    // Create the mapping.
    field_op.create_copy_mapping();

    // Do the transfer.
    field_op.copy();

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

TEUCHOS_UNIT_TEST( CopyOperator, Separate_Split_Transfer_Test )
{
    // Setup some communicators.
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;
    RCP_Communicator global_comm = getDefaultComm<int>();
    std::vector<int> source_ranks;
    std::vector<int> target_ranks;
    for ( int n = 0; n < global_comm->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    source_ranks.push_back(n);
	}
	else
	{
	    target_ranks.push_back(n);
	}
    }
    Teuchos::ArrayView<int> source_ranks_view( source_ranks );
    Teuchos::ArrayView<int> target_ranks_view( target_ranks );

    RCP_Communicator source_comm = 
	global_comm->createSubcommunicator( source_ranks_view );
    RCP_Communicator target_comm = 
	global_comm->createSubcommunicator( target_ranks_view );

    // Declare the target containers.
    Teuchos::RCP<Data_Container> source_container;
    Teuchos::RCP<Data_Container> target_container;

    // Declare the interfaces implementations on the global communicator.
    Teuchos::RCP<DataSource<double,int,double,3> > tds;
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt;

    // Create a data container instance for checking the data under the source
    // interface and create an instance of the source interface.
    if ( source_comm != Teuchos::null )
    {
	source_container = Teuchos::rcp(new Data_Container);

	tds = Teuchos::rcp(
	    new test_DataSource<double,int,double,3>(source_container,
						     source_comm));
    }

    // Create a data container instance for checking the data under the target
    // interface and create an instance of the target interface.
    if ( target_comm != Teuchos::null )
    {
	target_container = Teuchos::rcp(new Data_Container);

	tdt = Teuchos::rcp(
	    new test_DataTarget<double,int,double,3>(target_container,
						     target_comm));
    }

    // Create a distributed field for these interfaces to be transferred.
    global_comm->barrier();
    CopyOperator<double,int,double,3> field_op( global_comm,
						"DISTRIBUTED_TEST_FIELD", 
						"DISTRIBUTED_TEST_FIELD",
						tds, 
						tdt );

    // Create the mapping.
    field_op.create_copy_mapping();

    // Do the transfer.
    field_op.copy();

    // Check the transferred data under the target interface.
    if ( target_comm != Teuchos::null )
    {
	TEST_ASSERT( target_container->get_distributed_data().size() == 5 );
	int myRank = target_comm->getRank();
	int mySize = target_comm->getSize();
	int flippedRank = mySize-myRank-1;
	for (int i = 0; i < 5; ++i)
	{
	    TEST_ASSERT( target_container->get_distributed_data()[i]
			 == 1.0*flippedRank );
	}
    }
}

TEUCHOS_UNIT_TEST( CopyOperator, Overlap_Split_Transfer_Test )
{
    // Setup some communicators.
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;
    RCP_Communicator global_comm = getDefaultComm<int>();
    std::vector<int> source_ranks;
    std::vector<int> target_ranks;
    for ( int n = 0; n < global_comm->getSize(); ++n )
    {
	if ( n % 2 == 0 || n % 3 == 0 )
	{
	    source_ranks.push_back(n);
	}
	if ( n % 2 == 1 )
	{
	    target_ranks.push_back(n);
	}
    }
    Teuchos::ArrayView<int> source_ranks_view( source_ranks );
    Teuchos::ArrayView<int> target_ranks_view( target_ranks );

    RCP_Communicator source_comm = 
	global_comm->createSubcommunicator( source_ranks_view );
    RCP_Communicator target_comm = 
	global_comm->createSubcommunicator( target_ranks_view );

    // Declare the target containers.
    Teuchos::RCP<Data_Container> source_container;
    Teuchos::RCP<Data_Container> target_container;

    // Declare the interfaces implementations on the global communicator.
    Teuchos::RCP<DataSource<double,int,double,3> > tds;
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt;

    // Create a data container instance for checking the data under the source
    // interface and create an instance of the source interface.
    if ( source_comm != Teuchos::null )
    {
	source_container = Teuchos::rcp(new Data_Container);

	tds = Teuchos::rcp(
	    new test_DataSource<double,int,double,3>(source_container,
						     source_comm));
    }

    // Create a data container instance for checking the data under the target
    // interface and create an instance of the target interface.
    if ( target_comm != Teuchos::null )
    {
	target_container = Teuchos::rcp(new Data_Container);

	tdt = Teuchos::rcp(
	    new test_DataTarget<double,int,double,3>(target_container,
						     target_comm));
    }

    // Create a distributed field for these interfaces to be transferred.
    global_comm->barrier();
    CopyOperator<double,int,double,3> field_op( global_comm,
						"DISTRIBUTED_TEST_FIELD", 
						"DISTRIBUTED_TEST_FIELD",
						tds, 
						tdt );

    // Create the mapping.
    field_op.create_copy_mapping();

    // Do the transfer.
    field_op.copy();

    // Check the transferred data under the target interface.
    if ( target_comm != Teuchos::null )
    {
	TEST_ASSERT( target_container->get_distributed_data().size() == 5 );
	int myRank = target_comm->getRank();
	int mySize = target_comm->getSize();
	int flippedRank = mySize-myRank-1;
	for (int i = 0; i < 5; ++i)
	{
	    TEST_ASSERT( target_container->get_distributed_data()[i]
			 == 1.0*flippedRank );
	}
    }
}

TEUCHOS_UNIT_TEST( CopyOperator, Source_Subcomm_Transfer_Test )
{
    // Setup some communicators.
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;
    RCP_Communicator global_comm = getDefaultComm<int>();
    std::vector<int> source_ranks;
    for ( int n = 0; n < global_comm->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    source_ranks.push_back(n);
	}
    }
    Teuchos::ArrayView<int> source_ranks_view( source_ranks );

    RCP_Communicator source_comm = 
	global_comm->createSubcommunicator( source_ranks_view );
    RCP_Communicator target_comm = global_comm;

    // Declare the target containers.
    Teuchos::RCP<Data_Container> source_container;
    Teuchos::RCP<Data_Container> target_container;

    // Declare the interfaces implementations on the global communicator.
    Teuchos::RCP<DataSource<double,int,double,3> > tds;
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt;

    // Create a data container instance for checking the data under the source
    // interface and create an instance of the source interface.
    if ( source_comm != Teuchos::null )
    {
	source_container = Teuchos::rcp(new Data_Container);

	tds = Teuchos::rcp(
	    new test_DataSource<double,int,double,3>(source_container,
						     source_comm));
    }

    // Create a data container instance for checking the data under the target
    // interface and create an instance of the target interface.
    if ( target_comm != Teuchos::null )
    {
	target_container = Teuchos::rcp(new Data_Container);

	tdt = Teuchos::rcp(
	    new test_DataTarget<double,int,double,3>(target_container,
						     target_comm));
    }

    // Create a distributed field for these interfaces to be transferred.
    global_comm->barrier();
    CopyOperator<double,int,double,3> field_op( global_comm,
						"DISTRIBUTED_TEST_FIELD", 
						"DISTRIBUTED_TEST_FIELD",
						tds, 
						tdt );

    // Create the mapping.
    field_op.create_copy_mapping();

    // Do the transfer.
    field_op.copy();

    // Check the transferred data under the target interface.
    if ( target_comm != Teuchos::null )
    {
	if ( global_comm->getRank() % 2 == 0 )
	{
	    TEST_ASSERT( target_container->get_distributed_data().size() == 5 );
	    int myRank = source_comm->getRank();
	    int mySize = source_comm->getSize();
	    int flippedRank = mySize-myRank-1;
	    for (int i = 0; i < 5; ++i)
	    {
		TEST_ASSERT( target_container->get_distributed_data()[i]
			     == 1.0*flippedRank );
	    }
	}
	else
	{
	    TEST_ASSERT( target_container->get_distributed_data().size() == 5 );
	    for (int i = 0; i < 5; ++i)
	    {
		TEST_ASSERT( target_container->get_distributed_data()[i] == 0.0 );
	    }
	}
    }
}

TEUCHOS_UNIT_TEST( CopyOperator, Target_Subcomm_Transfer_Test )
{
    // Setup some communicators.
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;
    RCP_Communicator global_comm = getDefaultComm<int>();
    std::vector<int> target_ranks;
    for ( int n = 0; n < global_comm->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    target_ranks.push_back(n);
	}
    }
    Teuchos::ArrayView<int> target_ranks_view( target_ranks );

    RCP_Communicator source_comm = global_comm;
    RCP_Communicator target_comm = 
	global_comm->createSubcommunicator( target_ranks_view );


    // Declare the target containers.
    Teuchos::RCP<Data_Container> source_container;
    Teuchos::RCP<Data_Container> target_container;

    // Declare the interfaces implementations on the global communicator.
    Teuchos::RCP<DataSource<double,int,double,3> > tds;
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt;

    // Create a data container instance for checking the data under the source
    // interface and create an instance of the source interface.
    if ( source_comm != Teuchos::null )
    {
	source_container = Teuchos::rcp(new Data_Container);

	tds = Teuchos::rcp(
	    new test_DataSource<double,int,double,3>(source_container,
						     source_comm));
    }

    // Create a data container instance for checking the data under the target
    // interface and create an instance of the target interface.
    if ( target_comm != Teuchos::null )
    {
	target_container = Teuchos::rcp(new Data_Container);

	tdt = Teuchos::rcp(
	    new test_DataTarget<double,int,double,3>(target_container,
						     target_comm));
    }

    // Create a distributed field for these interfaces to be transferred.
    global_comm->barrier();
    CopyOperator<double,int,double,3> field_op( global_comm,
						"DISTRIBUTED_TEST_FIELD", 
						"DISTRIBUTED_TEST_FIELD",
						tds, 
						tdt );

    // Create the mapping.
    field_op.create_copy_mapping();

    // Do the transfer.
    field_op.copy();

    // Check the transferred data under the target interface.
    if ( target_comm != Teuchos::null )
    {
	if ( global_comm->getRank() % 2 == 0 )
	{
	    TEST_ASSERT( target_container->get_distributed_data().size() == 5 );
	    int myRank = source_comm->getRank();
	    int mySize = source_comm->getSize();
	    int flippedRank = mySize-myRank-1;
	    for (int i = 0; i < 5; ++i)
	    {
		TEST_ASSERT( target_container->get_distributed_data()[i]
			     == 1.0*flippedRank );
	    }
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace Coupler

//---------------------------------------------------------------------------//
//                        end of tstCopyOperator.cpp
//---------------------------------------------------------------------------//
