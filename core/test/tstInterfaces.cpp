//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstInterfaces.cpp
 * \author Stuart Slattery
 * \date   Thu Dec 01 16:50:04 2011
 * \brief  Unit tests for the data transfer pure virtual interfaces.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DataTransferKit_Point.hpp>
#include <DataTransferKit_DataSource.hpp>
#include <DataTransferKit_DataTarget.hpp>

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"
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

    typedef DataTransferKit::Point<int,double>     PointType;

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
	return distributed_points;
    }

    Teuchos::ArrayView<double> get_distributed_data()
    {
	distributed_data.resize(1);
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

namespace DataTransferKit {

// transfer data source implementation - this implementation specifies double
// as the data type
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class test_DataSource 
    : public DataSource<DataType_T, HandleType_T, CoordinateType_T>
{
  private:

    std::vector<double> private_data;

  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef Point<HandleType,CoordinateType>         PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;

    test_DataSource()
    { /* ... */ }

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

    bool is_local_point(const PointType &point)
    {
	bool return_val = false;

	if ( point.x() > 0 && point.y() > 0 && point.z() > 0 )
	{
	    return_val = true;
	}

	return return_val;
    }

    const Teuchos::ArrayView<double> 
    get_source_data(const std::string &field_name)
    {
	Teuchos::ArrayView<double> return_view;

	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    private_data.clear();
	    private_data.resize(1);
	    private_data[0] = 1.0*getDefaultComm<OrdinalType>()->getRank();
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
	    return_val = 1.0*getDefaultComm<OrdinalType>()->getRank();
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
    std::vector<PointType> points;

  public:

    test_DataTarget(Teuchos::RCP<Data_Container> _container)
	: container(_container)
    { /* ... */ }

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
	    PointType local_point(1, 1.0, 1.0, 1.0);
	    std::vector<PointType> local_points(1, local_point);
	    points = local_points;
	    return_view = Teuchos::ArrayView<PointType>(points);
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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace DataTransferKit {

TEUCHOS_UNIT_TEST( Transfer_DataSource, source_interface_test )
{
    typedef Point<int,double>     PointType;

    // create an instance of the source interface.
    Teuchos::RCP<DataSource<double,int,double> > source_iface = 
	Teuchos::rcp(new test_DataSource<double,int,double>);

    // test the interface methods.
    TEST_ASSERT( source_iface->get_source_comm()->getRank() 
		 == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( source_iface->get_source_comm()->getSize() 
		 == getDefaultComm<int>()->getSize() );

    TEST_ASSERT( source_iface->is_field_supported("DISTRIBUTED_TEST_FIELD") );
    TEST_ASSERT( source_iface->is_field_supported("SCALAR_TEST_FIELD") );
    TEST_ASSERT( !source_iface->is_field_supported("FOO_TEST_FIELD") );

    PointType positive_point(1, 1.0, 1.0, 1.0);
    PointType negative_point(1, -1.0, -1.0, -1.0);
    TEST_ASSERT( source_iface->is_local_point(positive_point) );
    TEST_ASSERT( !source_iface->is_local_point(negative_point) );

    Teuchos::ArrayView<double> data_view;
    data_view = source_iface->get_source_data("FOO_TEST_FIELD");
    TEST_ASSERT( data_view.is_null() );
    data_view = source_iface->get_source_data("DISTRIBUTED_TEST_FIELD");
    TEST_ASSERT( data_view.size() == 1);
    TEST_ASSERT( data_view[0] == 1.0*getDefaultComm<int>()->getRank() );

    TEST_ASSERT( source_iface->get_global_source_data("FOO_TEST_FIELD") == 0.0 );
    TEST_ASSERT( source_iface->get_global_source_data("SCALAR_TEST_FIELD") == 
		 1.0*getDefaultComm<int>()->getRank() );
}

TEUCHOS_UNIT_TEST( Transfer_DataTarget, target_interface_test )
{
    typedef Point<int,double>     PointType;

    // create a data container instance.
    Teuchos::RCP<Data_Container> container = Teuchos::rcp(new Data_Container);

    // create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double> > target_iface = 
	Teuchos::rcp(
	    new test_DataTarget<double,int,double>(container));

    // test the interface methods.
    TEST_ASSERT( target_iface->get_target_comm()->getRank() 
		 == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( target_iface->get_target_comm()->getSize() 
		 == getDefaultComm<int>()->getSize() );

    TEST_ASSERT( target_iface->is_field_supported("DISTRIBUTED_TEST_FIELD") );
    TEST_ASSERT( target_iface->is_field_supported("SCALAR_TEST_FIELD") );
    TEST_ASSERT( !target_iface->is_field_supported("FOO_TEST_FIELD") );

    Teuchos::ArrayView<PointType> points_view;
    points_view = target_iface->get_target_points("FOO_TEST_FIELD");
    TEST_ASSERT( points_view.size() == 0 );
    points_view = target_iface->get_target_points("DISTRIBUTED_TEST_FIELD");
    TEST_ASSERT( points_view.size() == 1 );
    TEST_ASSERT( points_view[0].handle() == 1 );
    TEST_ASSERT( points_view[0].x() == 1.0 );
    TEST_ASSERT( points_view[0].y() == 1.0 );
    TEST_ASSERT( points_view[0].z() == 1.0 );

    Teuchos::ArrayView<double> data_view;
    data_view = target_iface->get_target_data_space("FOO_TEST_FIELD");
    TEST_ASSERT( data_view.is_null() );
    data_view = target_iface->get_target_data_space("DISTRIBUTED_TEST_FIELD");
    data_view[0] = 1.0*getDefaultComm<int>()->getRank();
    TEST_ASSERT( container->get_distributed_data().size() == 1 );
    TEST_ASSERT( container->get_distributed_data()[0] == 
		 1.0*getDefaultComm<int>()->getRank() );

    double global_scalar = 1.0*getDefaultComm<int>()->getRank();
    target_iface->set_global_target_data("FOO_TEST_FIELD", global_scalar);
    target_iface->set_global_target_data("SCALAR_TEST_FIELD", global_scalar);
    TEST_ASSERT( container->get_scalar_data() == 
		 1.0*getDefaultComm<int>()->getRank() );
}

TEUCHOS_UNIT_TEST( Transfer_DataSource, simple_coupling_test )
{
    typedef Point<int,double>     PointType;

    // create an instance of the source interface.
    Teuchos::RCP<DataSource<double,int,double> > source_iface = 
	Teuchos::rcp(new test_DataSource<double,int,double>);

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> container = Teuchos::rcp(new Data_Container);

    // create an instance of the target interface.
    Teuchos::RCP<DataTarget<double,int,double> > target_iface = 
	Teuchos::rcp( new test_DataTarget<double,int,double>(container));

    // Check that the field is supported.
    TEST_ASSERT( source_iface->is_field_supported("DISTRIBUTED_TEST_FIELD") &&
		 target_iface->is_field_supported("DISTRIBUTED_TEST_FIELD") );

    // Test a target point in the source
    TEST_ASSERT( source_iface->is_local_point( 
		     target_iface->get_target_points(
			 "DISTRIBUTED_TEST_FIELD")[0] ) );

    // Transfer data from the source to the target.
    Teuchos::ArrayView<double> source_view = 
	source_iface->get_source_data("DISTRIBUTED_TEST_FIELD");

    Teuchos::ArrayView<double> target_view = 
	target_iface->get_target_data_space("DISTRIBUTED_TEST_FIELD");

    target_view.assign(source_view);

    TEST_ASSERT( container->get_distributed_data().size() == 1 );
    TEST_ASSERT( container->get_distributed_data()[0] == 
		 1.0*getDefaultComm<int>()->getRank() );
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
//                        end of tstInterfaces.cpp
//---------------------------------------------------------------------------//
