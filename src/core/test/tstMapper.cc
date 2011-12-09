//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstMapper.cc
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:31:19 2011
 * \brief  Unit tests for the mapper operator.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "../Transfer_Data_Source.hh"
#include "../Transfer_Data_Target.hh"
#include "../Transfer_Data_Field.hh"
#include "../Mapper_Tpetra.hh"

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
    Teuchos::ArrayView<PointType > distributed_points;

  public:

    Data_Container()
    { /* ... */ }

    ~Data_Container()
    { /* ... */ }

    void set_distributed_points(Teuchos::ArrayView<PointType > points)
    {
	distributed_points = points;
    }

    Teuchos::ArrayView<PointType > get_distributed_points()
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

	OrdinalType my_rank = getDefaultComm<OrdinalType>()->getRank();

	if ( point.x() == 1.0*my_rank &&
	     point.y() == 2.0*my_rank &&
	     point.z() == 3.0*my_rank )
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

    std::vector<PointType> points;
    Teuchos::RCP<Data_Container> container;

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
	    OrdinalType my_rank = getDefaultComm<OrdinalType>()->getRank();
	    PointType local_point(my_rank, 
				  1.0*my_rank, 
				  2.0*my_rank, 
				  3.0*my_rank);
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

TEUCHOS_UNIT_TEST( Mapper, mirror_test )
{
    // Create a data container instance.
    Teuchos::RCP<Data_Container> container = Teuchos::rcp(new Data_Container);

    // Create an instance of the source interface.
    Teuchos::RCP<Transfer_Data_Source<double,int,double> > tds = 
	Teuchos::rcp(new test_Transfer_Data_Source<double,int,double>());

    // Create an instance of the target interface.
    Teuchos::RCP<Transfer_Data_Target<double,int,double> > tdt = 
	Teuchos::rcp(new test_Transfer_Data_Target<double,int,double>(container));

    // Create a distributed field for these interfaces to be transferred.
    Teuchos::RCP<Transfer_Data_Field<double,int,double> > field = Teuchos::rcp(
	new Transfer_Data_Field<double,int,double>("DISTRIBUTED_TEST_FIELD", tds, tdt));

    // Create a mapper and generate the map.
    Mapper<double,int,double> mapper;
    mapper.map(getDefaultComm<int>(), field);

    // Check the contents of the map - everyone should be sending and
    // receiving only from themselves.

}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                        end of tstMapper.cc
//---------------------------------------------------------------------------//
