//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstInterfaces.cc
 * \author Stuart Slattery
 * \date   Thu Dec 01 16:50:04 2011
 * \brief  Unit tests for the data transfer pure virtual interfaces.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Mesh_Point.hpp>
#include "../Transfer_Data_Source.hh"
#include "../Transfer_Data_Target.hh"

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
#ifdef COMM_MPI
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

    typedef mesh::Point<int,double>     PointType;

  private:

    double scalar_data;
    Teuchos::ArrayView<double> distributed_data;
    Teuchos::ArrayView<PointType> distributed_points;

  public:

    Data_Container()
    { /* ... */ }

    ~Data_Container()
    { /* ... */ }

    void set_distributed_points(Teuchos::ArrayView<PointType> points)
    {
	distributed_points = points;
    }

    Teuchos::ArrayView<PointType> get_distributed_points()
    {
	return distributed_points;
    }

    void set_distributed_data(Teuchos::ArrayView<double> data)
    {
	distributed_data = data;
    }

    Teuchos::ArrayView<double> get_distributed_data()
    {
	return distributed_data;
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

namespace coupler {

// transfer data source implementation - this implementation specifies double
// as the data type
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class test_Transfer_Data_Source 
    : public Transfer_Data_Source<DataType_T, HandleType_T, CoordinateType_T>
{
  private:

    std::vector<double> private_data;

  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef mesh::Point<HandleType,CoordinateType>   PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;

    test_Transfer_Data_Source()
    { /* ... */ }

    ~test_Transfer_Data_Source()
    { /* ... */ }

    const RCP_Communicator comm()
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

	if ( point.x() > 0 && point.y() > 0 && point.z() > 0 )
	{
	    return_val = true;
	}

	return return_val;
    }

    const Teuchos::ArrayView<double> send_data(const std::string &field_name)
    {
	Teuchos::ArrayView<double> return_view;

	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    std::vector<double> local_data(1, 1.0);
	    private_data = local_data;
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
	    return_val = 1.0;
	}

	return return_val;
    }
};

//---------------------------------------------------------------------------//

// transfer data target implementation - this implementation specifies double
// as the data type
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class test_Transfer_Data_Target 
    : public Transfer_Data_Target<DataType_T, HandleType_T, CoordinateType_T>
{
  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef mesh::Point<HandleType,CoordinateType>   PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;

  private:

    Teuchos::RCP<Data_Container> container;
    std::vector<PointType> points;

  public:

    test_Transfer_Data_Target(Teuchos::RCP<Data_Container> _container)
	: container(_container)
    { /* ... */ }

    ~test_Transfer_Data_Target()
    { /* ... */ }

    const RCP_Communicator comm()
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
	    PointType local_point(1, 1.0, 1.0, 1.0);
	    std::vector<PointType> local_points(1, local_point);
	    points = local_points;
	    return_view = Teuchos::ArrayView<PointType>(points);
	}

	return return_view;
    }

    void receive_data(const std::string &field_name,
		      const Teuchos::ArrayView<DataType> &data)
    {
	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    container->set_distributed_data(data);
	}
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

} // end namespace coupler

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace coupler {

TEUCHOS_UNIT_TEST( Transfer_Data_Source, source_interface_test )
{
    typedef mesh::Point<int,double>     PointType;

    // create an instance of the source interface.
    Teuchos::RCP<Transfer_Data_Source<double,int,double> > source_iface = 
	Teuchos::rcp(new test_Transfer_Data_Source<double,int,double>);

    // test the interface methods.
    TEST_ASSERT( source_iface->comm()->getSize() > 
		 source_iface->comm()->getRank() );
    std::cout << "Source comm size: " 
	      << source_iface->comm()->getSize() << std::endl;
    std::cout << "On rank: " << source_iface->comm()->getRank() << std::endl;

    TEST_ASSERT( source_iface->field_supported("DISTRIBUTED_TEST_FIELD") );
    TEST_ASSERT( source_iface->field_supported("SCALAR_TEST_FIELD") );
    TEST_ASSERT( !source_iface->field_supported("FOO_TEST_FIELD") );

    PointType positive_point(1, 1.0, 1.0, 1.0);
    PointType negative_point(1, -1.0, -1.0, -1.0);
    TEST_ASSERT( source_iface->get_points(positive_point) );
    TEST_ASSERT( !source_iface->get_points(negative_point) );

    Teuchos::ArrayView<double> data_view;
    data_view = source_iface->send_data("FOO_TEST_FIELD");
    TEST_ASSERT( data_view.size() == 0 );
    data_view = source_iface->send_data("DISTRIBUTED_TEST_FIELD");
    TEST_ASSERT( data_view.size() == 1);
    TEST_ASSERT( data_view[0] == 1.0 );

    TEST_ASSERT( source_iface->set_global_data("FOO_TEST_FIELD") == 0.0 );
    TEST_ASSERT( source_iface->set_global_data("SCALAR_TEST_FIELD") == 1.0 );
}

TEUCHOS_UNIT_TEST( Transfer_Data_Target, target_interface_test )
{
    typedef mesh::Point<int,double>     PointType;

    // create a data container instance.
    Teuchos::RCP<Data_Container> container = Teuchos::rcp(new Data_Container);

    // create an instance of the target interface.
    Teuchos::RCP<Transfer_Data_Target<double,int,double> > target_iface = 
	Teuchos::rcp(
	    new test_Transfer_Data_Target<double,int,double>(container));

    // test the interface methods.
    TEST_ASSERT( target_iface->comm()->getSize() > 
		 target_iface->comm()->getRank() );
    std::cout << "Target comm size: " << target_iface->comm()->getSize() << std::endl;
    std::cout << "On rank: " << target_iface->comm()->getRank() << std::endl;

    TEST_ASSERT( target_iface->field_supported("DISTRIBUTED_TEST_FIELD") );
    TEST_ASSERT( target_iface->field_supported("SCALAR_TEST_FIELD") );
    TEST_ASSERT( !target_iface->field_supported("FOO_TEST_FIELD") );

    Teuchos::ArrayView<PointType> points_view;
    points_view = target_iface->set_points("FOO_TEST_FIELD");
    TEST_ASSERT( points_view.size() == 0 );
    points_view = target_iface->set_points("DISTRIBUTED_TEST_FIELD");
    TEST_ASSERT( points_view.size() == 1 );
    TEST_ASSERT( points_view[0].handle() == 1 );
    TEST_ASSERT( points_view[0].x() == 1.0 );
    TEST_ASSERT( points_view[0].y() == 1.0 );
    TEST_ASSERT( points_view[0].z() == 1.0 );

    std::vector<double> data_to_receive(1, 1.0);
    Teuchos::ArrayView<double> data_view(data_to_receive);
    target_iface->receive_data("FOO_TEST_FIELD", data_view);
    TEST_ASSERT( container->get_distributed_data().size() == 0 );
    TEST_ASSERT( container->get_distributed_points().size() == 0 );
    target_iface->receive_data("DISTRIBUTED_TEST_FIELD", data_view);
    TEST_ASSERT( container->get_distributed_data().size() == 1 );
    TEST_ASSERT( container->get_distributed_data()[0] == 1.0 );

    double global_scalar = 1.0;
    target_iface->get_global_data("FOO_TEST_FIELD", global_scalar);
    TEST_ASSERT( container->get_scalar_data() != 1.0 );
    target_iface->get_global_data("SCALAR_TEST_FIELD", global_scalar);
    TEST_ASSERT( container->get_scalar_data() == 1.0 );
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                        end of tstInterfaces.cc
//---------------------------------------------------------------------------//
