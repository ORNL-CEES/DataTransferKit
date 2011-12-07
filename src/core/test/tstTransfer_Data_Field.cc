//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstTransfer_Data_Field.cc
 * \author Stuart Slattery
 * \date   Fri Nov 18 14:43:10 2011
 * \brief  Transfer_Data_Field unit tests
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "comm/global.hh"

#include <Mesh_Point.hpp>

#include "../Transfer_Map.hh"
#include "../Transfer_Data_Source.hh"
#include "../Transfer_Data_Target.hh"
#include "../Transfer_Data_Field.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_UnitTestHarness.hpp"

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
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef mesh::Point<HandleType,CoordinateType>   PointType;

    test_Transfer_Data_Source()
    { /* ... */ }

    ~test_Transfer_Data_Source()
    { /* ... */ }

    void register_comm(Communicator &comm)
    {
#ifdef COMM_MPI
	comm = MPI_COMM_WORLD;
#else
	comm = 1;
#endif
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

    bool get_points(PointType &point)
    {
	bool return_val = false;

	if ( point.x() > 0 && point.y() > 0 && point.z() > 0 )
	{
	    return_val = true;
	}

	return return_val;
    }

    Teuchos::ArrayView<double> send_data(
	const std::string &field_name,
	const Teuchos::ArrayView<PointType> &points)
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
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef mesh::Point<HandleType,CoordinateType>   PointType;

  private:

    std::vector<PointType> points;

  public:

    test_Transfer_Data_Target()
    { /* ... */ }

    ~test_Transfer_Data_Target()
    { /* ... */ }

    void register_comm(Communicator &comm)
    {
#ifdef COMM_MPI
	comm = MPI_COMM_WORLD;
#else
	comm = 1;
#endif
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

    Teuchos::ArrayView<PointType> set_points(const std::string &field_name)
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
		      const Teuchos::ArrayView<PointType> &points,
		      const Teuchos::ArrayView<DataType> &data)
    { /* ... */ }

    void get_global_data(const std::string &field_name,
			 const DataType &data)
    { /* ... */ }
};

} // end namespace coupler

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace coupler {

TEUCHOS_UNIT_TEST( Transfer_Data_Field, distributed_field_test )
{
    // Create an instance of the source interface.
    Teuchos::RCP<Transfer_Data_Source<double,int,double> > tds = 
	Teuchos::rcp(new test_Transfer_Data_Source<double,int,double>());

    // Create an instance of the target interface.
    Teuchos::RCP<Transfer_Data_Target<double,int,double> > tdt = 
	Teuchos::rcp(new test_Transfer_Data_Target<double,int,double>());

    // Create a distributed field for these interfaces to be transferred.
    Transfer_Data_Field<double,int,double> field("DISTRIBUTED_TEST_FIELD", tds, tdt);

    // Add a transfer map to the field.
    TEST_ASSERT( !field.is_mapped() );
    Teuchos::RCP<Transfer_Map> map = Teuchos::rcp(new Transfer_Map());
    field.set_map(map);

    // Test the functionality.
    TEST_ASSERT( field.name() == "DISTRIBUTED_TEST_FIELD" );
    TEST_ASSERT( field.source() == tds );
    TEST_ASSERT( field.target() == tdt );
    TEST_ASSERT( field.get_map() == map );
    TEST_ASSERT( !field.is_scalar() );
    TEST_ASSERT( field.is_mapped() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Transfer_Data_Field, scalar_field_test )
{
    // Create an instance of the source interface.
    Teuchos::RCP<Transfer_Data_Source<double,int,double> > tds = 
	Teuchos::rcp(new test_Transfer_Data_Source<double,int,double>());

    // Create an instance of the target interface.
    Teuchos::RCP<Transfer_Data_Target<double,int,double> > tdt = 
	Teuchos::rcp(new test_Transfer_Data_Target<double,int,double>());

    // Create a scalar field for these interfaces to be transferred.
    Transfer_Data_Field<double,int,double> field("SCALAR_TEST_FIELD", tds, tdt, true);

    // Test the functionality.
    TEST_ASSERT( field.name() == "SCALAR_TEST_FIELD" );
    TEST_ASSERT( field.source() == tds );
    TEST_ASSERT( field.target() == tdt );
    TEST_ASSERT( field.is_scalar() );
    TEST_ASSERT( !field.is_mapped() );
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                        end of tstTransfer_Data_Field.cc
//---------------------------------------------------------------------------//
