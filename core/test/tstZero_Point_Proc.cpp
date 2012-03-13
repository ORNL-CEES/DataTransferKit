//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstZero_Point_Proc.cpp
 * \author Stuart Slattery
 * \date   Fri Nov 18 14:43:10 2011
 * \brief  CopyOperator unit tests
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "Coupler_Point.hpp"
#include "Coupler_DataSource.hpp"
#include "Coupler_DataTarget.hpp"
#include "Coupler_CopyOperator.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

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

  public:

    test_DataSource(Teuchos::RCP<Data_Container> _container)
	: container(_container)
    { 
	myRank = getDefaultComm<OrdinalType>()->getRank();
	mySize = getDefaultComm<OrdinalType>()->getSize();
    }

    ~test_DataSource()
    { /* ... */ }

    RCP_Communicator get_source_comm()
    {
	return getDefaultComm<OrdinalType>();
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

	// We only want to find the first 4 of the 5 points from the simpler
	// test. This will give maps and vectors of different sizes.
	if ( (int) local_points.size() != 4 &&
	     test_point.getCoords()[0] == 1.0*(mySize-myRank-1) &&
	     test_point.getCoords()[1] == 2.0*(mySize-myRank-1) &&
	     test_point.getCoords()[2] == 3.0*(mySize-myRank-1) )
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

  public:

    test_DataTarget(Teuchos::RCP<Data_Container> _container)
	: container(_container)
    { 
	myRank = getDefaultComm<OrdinalType>()->getRank();
	mySize = getDefaultComm<OrdinalType>()->getSize();
    }

    ~test_DataTarget()
    { /* ... */ }

    RCP_Communicator get_target_comm()
    {
	return getDefaultComm<OrdinalType>();
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
	    // Make proc 0 a zero point proc.
	    if ( myRank != 0 )
	    {
		for (int i = 0; i < 5; ++i)
		{
		    int new_handle = 5*myRank + i;
		    PointType new_point = point(new_handle, 1.0*myRank, 
						2.0*myRank, 3.0*myRank);
		    local_points.push_back(new_point);
		}
		return_view = Teuchos::ArrayView<PointType>(local_points);
	    }
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
	Teuchos::rcp(new test_DataSource<double,int,double,3>(source_container));

    // Create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double,3> > tdt = 
	Teuchos::rcp(new test_DataTarget<double,int,double,3>(target_container));

    // Create a distributed field for these interfaces to be transferred.
    CopyOperator<double,int,double,3> field_op( getDefaultComm<int>(),
						"DISTRIBUTED_TEST_FIELD",
						"DISTRIBUTED_TEST_FIELD", 
						tds, 
						tdt );

    // Create the mapping.
    field_op.create_copy_mapping();

    // Do the transfer.
    field_op.copy();

    // Check the points under the source interface to the verify communication
    // pattern of sending target points to the source.
    int myRank = getDefaultComm<int>()->getRank();
    int mySize = getDefaultComm<int>()->getSize();
    int flippedRank = mySize-myRank-1;

    TEST_ASSERT( field_op.is_mapped() );
    if ( flippedRank != 0 )
    {
	TEST_ASSERT( source_container->get_distributed_points().size() == 4 );
	for (int i = 0; i < 4; ++i)
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
    }

    // Check the Tpetra maps for both the source and the target.
    TEST_ASSERT( (int) field_op.source_map()->getGlobalNumElements() 
		 == (mySize-1)*4 );
    if ( flippedRank != 0 )
    {
	TEST_ASSERT( (int) field_op.source_map()->getNodeNumElements() == 4 );
	for (int i = 0; i < 4; ++i)
	{
	    TEST_ASSERT( field_op.source_map()->getNodeElementList()[i] 
			 == 5*flippedRank+i );
	}
    }

    TEST_ASSERT( (int) field_op.target_map()->getGlobalNumElements() 
		 == (mySize-1)*5 );
    if ( myRank != 0 )
    {
	TEST_ASSERT( (int) field_op.target_map()->getNodeNumElements() == 5 );
	for (int i = 0; i < 5; ++i)
	{
	    TEST_ASSERT( field_op.target_map()->getNodeElementList()[i] 
			 == 5*myRank+i );
	}

	// Check the transferred data under the target interface.
	TEST_ASSERT( target_container->get_distributed_data().size() == 5 );
	for (int i = 0; i < 4; ++i)
	{
	    TEST_ASSERT( target_container->get_distributed_data()[i]
			 == 1.0*flippedRank );
	}
    
	// The last element in each view should be zero because we told the
	// source to skip those points.
	TEST_ASSERT ( target_container->get_distributed_data()[4] == 0.0 );
    }
}

//---------------------------------------------------------------------------//

} // end namespace Coupler

//---------------------------------------------------------------------------//
//                        end of tstZero_Point_Proc.cpp
//---------------------------------------------------------------------------//
