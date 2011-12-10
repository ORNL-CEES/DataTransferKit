//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstTransfer_Data_Field.cpp
 * \author Stuart Slattery
 * \date   Fri Nov 18 14:43:10 2011
 * \brief  Transfer_Data_Field unit tests
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Mesh_Point.hpp>

#include <Coupler_Data_Source.hpp>
#include <Coupler_Data_Target.hpp>
#include <Coupler_Data_Field.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Tpetra_Map.hpp"

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

namespace coupler {

// transfer data source implementation - this implementation specifies double
// as the data type
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class test_Data_Source 
    : public Data_Source<DataType_T, HandleType_T, CoordinateType_T>
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

    test_Data_Source()
    { /* ... */ }

    ~test_Data_Source()
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
	    private_data.clear();
	    private_data.resize(1);
	    private_data[0] = 1.0*getDefaultComm<OrdinalType>()->getRank();
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
	    return_val = 1.0*getDefaultComm<OrdinalType>()->getRank();
	}

	return return_val;
    }
};

//---------------------------------------------------------------------------//
// transfer data target implementation - this implementation specifies double
// as the data type
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class test_Data_Target 
    : public Data_Target<DataType_T, HandleType_T, CoordinateType_T>
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

    test_Data_Target(Teuchos::RCP<Data_Container> _container)
	: container(_container)
    { /* ... */ }

    ~test_Data_Target()
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

} // end namespace coupler

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace coupler {

TEUCHOS_UNIT_TEST( Transfer_Data_Field, distributed_container_test )
{
    // Create an instance of the source interface.
    Teuchos::RCP<Data_Source<double,int,double> > tds = 
	Teuchos::rcp(new test_Data_Source<double,int,double>());

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> container = Teuchos::rcp(new Data_Container);

    // Create an instance of the target interface.
    Teuchos::RCP<Data_Target<double,int,double> > tdt = 
	Teuchos::rcp(new test_Data_Target<double,int,double>(container));

    // Create a distributed field for these interfaces to be transferred.
    Data_Field<double,int,double> field(getDefaultComm<int>(),
					"DISTRIBUTED_TEST_FIELD", 
					tds, 
					tdt);

    // Test the functionality.
    TEST_ASSERT( field.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( field.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( field.name() == "DISTRIBUTED_TEST_FIELD" );
    TEST_ASSERT( field.source() == tds );
    TEST_ASSERT( field.target() == tdt );
    TEST_ASSERT( !field.is_scalar() );
    TEST_ASSERT( field.is_mapped() );
}

TEUCHOS_UNIT_TEST( Transfer_Data_Field, scalar_container_test )
{
    // Create an instance of the source interface.
    Teuchos::RCP<Data_Source<double,int,double> > tds = 
	Teuchos::rcp(new test_Data_Source<double,int,double>());

    // create a data container instance for checking the data under the target
    // interface.
    Teuchos::RCP<Data_Container> container = Teuchos::rcp(new Data_Container);

    // Create an instance of the target interface.
    Teuchos::RCP<Data_Target<double,int,double> > tdt = 
	Teuchos::rcp(new test_Data_Target<double,int,double>(container));

    // Create a distributed field for these interfaces to be transferred.
    Data_Field<double,int,double> field(getDefaultComm<int>(),
					"DISTRIBUTED_TEST_FIELD", 
					tds, 
					tdt,
					true);

    // Test the functionality.
    TEST_ASSERT( field.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( field.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( field.name() == "DISTRIBUTED_TEST_FIELD" );
    TEST_ASSERT( field.source() == tds );
    TEST_ASSERT( field.target() == tdt );
    TEST_ASSERT( field.is_scalar() );
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                        end of tstTransfer_Data_Field.cpp
//---------------------------------------------------------------------------//
