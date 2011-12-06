//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstMessenger.cc
 * \author Stuart Slattery
 * \date   Thu Jun 02 09:10:58 2011
 * \brief  Tests the Messenger class.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template_c4_test.cc,v 1.7 2008/01/02 22:50:26 9te Exp $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "release/Release.hh"
#include "../Messenger.hh"
#include "../Transfer_Data_Source.hh"
#include "../Transfer_Data_Target.hh"
#include "../Transfer_Data_Field.hh"
#include "../Transfer_Map.hh"
#include "../Mapper.hh"
#include "Teuchos_RCP.hpp"

using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using coupler::Transfer_Data_Source;
using coupler::Transfer_Data_Target;
using coupler::Transfer_Data_Field;
using coupler::Transfer_Map;
using coupler::Mapper;
using coupler::Messenger;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
nemesis::Communicator_t get_comm_world()
{
#ifdef COMM_MPI
    return MPI_COMM_WORLD;
#else
    return 1;
#endif
}

//---------------------------------------------------------------------------//
// DATA CLASS
//---------------------------------------------------------------------------//

class Data_Container
{
  private:

    double scalar_data;
    std::vector<double> distributed_data;
    std::vector<int> distributed_handles;

  public:

    Data_Container()
    { /* ... */ }

    ~Data_Container()
    { /* ... */ }

    void set_distributed_handles(std::vector<int> handles)
    {
	distributed_handles = handles;
    }

    std::vector<int> get_distributed_handles()
    {
	return distributed_handles;
    }

    void set_distributed_data(std::vector<double> data)
    {
	distributed_data = data;
    }

    std::vector<double> get_distributed_data()
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

// transfer data source implementation - this implementation specifies double
// as the data type
template<class DataType_T>
class test_Transfer_Data_Source : public Transfer_Data_Source<DataType_T>
{
  public:

    //@{
    //! Useful typedefs.
    typedef double                                   DataType;
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    //@}

    /*!
     * \brief Constructor.
     */
    test_Transfer_Data_Source()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~test_Transfer_Data_Source()
    { /* ... */ }

    /*!
     * \brief Register communicator object.
     * \param comm The communicator for this physics.
     */
    void register_comm(Communicator &comm)
    {
#ifdef COMM_MPI
	comm = MPI_COMM_WORLD;
#else
	comm = 1;
#endif
    }

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported. 
     * \param field_name The name of the field for which support is being
     * checked.
     */
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

    /*! 
     * \brief Given (x,y,z) coordinates and an associated globally unique
     * handle, return true if the point is in the local domain, false if not.
     * \param handle The globally unique handle associated with the point.
     * \param x X coordinate.
     * \param y Y coordinate.
     * \param z Z coordinate.
     */
    bool get_points(HandleType handle,
		    CoordinateType x, 
		    CoordinateType y,
		    CoordinateType z)
    {
	bool return_val = false;

	if ( x == 1.0*nemesis::node() )
	{
	    return_val = true;
	}

	return return_val;
    }

    /*! 
     * \brief Given an entity handle, send the field data associated with that
     * handle. 
     * \param field_name The name of the field to send data from.
     * \param handles The enitity handles for the data being sent.
     * \param data The data being sent.
     */
    void send_data(const std::string &field_name,
		   const std::vector<HandleType> &handles,
		   std::vector<DataType> &data)
    {
	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    std::vector<double> local_data(1, 1.0*nemesis::node() );
	    data = local_data;
	}
    }

    /*!
     * \brief Given a field, set a global data element to be be sent to a
     * target.
     * \param field_name The name of the field to send data from.
     * \param data The global data element.
     */
    void set_global_data(const std::string &field_name,
			 DataType &data)
    {
	if ( field_name == "SCALAR_TEST_FIELD" )
	{
	    data = nemesis::node();
	}
    }
};

//---------------------------------------------------------------------------//

// transfer data target implementation - this implementation specifies double
// as the data type
template<class DataType_T>
class test_Transfer_Data_Target : public Transfer_Data_Target<DataType_T>
{
  private:

    Teuchos::RCP<Data_Container> container;

  public:

    //@{
    //! Useful typedefs.
    typedef double                                   DataType;
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    //@}

    /*!
     * \brief Constructor.
     */
    test_Transfer_Data_Target()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~test_Transfer_Data_Target()
    { /* ... */ }

    /*!
     * \brief Register communicator object.
     * \param comm The communicator for this physics.
     */
    void register_comm(Communicator &comm)
    {
#ifdef COMM_MPI
	comm = MPI_COMM_WORLD;
#else
	comm = 1;
#endif
    }

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported. 
     * \param field_name The name of the field for which support is being
     * checked.
     */
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

    /*!
     * \brief Set cartesian coordinates with a field. The coordinate
     * vector should be interleaved. The handle vector should consist of
     * globally unique handles. 
     * \param field_name The name of the field that the coordinates are being
     * registered with.
     * \param handles Point handle array.
     * \param coordinates Point coordinate array.
     */
    void set_points(const std::string &field_name,
		    std::vector<HandleType> &handles,
		    std::vector<CoordinateType> &coordinates)
    {
	if ( field_name = "DISTRIBUTED_TEST_FIELD" )
	{
	    std::vector<int> local_handles(1, nemesis::node() );
	    std::vector<double> local_coords(3, 1.0*nemesis::node() );

	    handles = local_handles;
	    coordinates = local_coordinates;
	}
    }

    /*! 
     * \brief Given an entity handle, receive the field data associated with
     * that handle. 
     * \param field_name The name of the field to receive data from.
     * \param handles The enitity handles for the data being received.
     * \param data The data being received.
     */
    void receive_data(const std::string &field_name,
		      const std::vector<HandleType> &handles,
		      const std::vector<DataType> &data)
    {
	if ( field_name = "DISTRIBUTED_TEST_FIELD" )
	{
	    container->set_distributed_handles(handles);
	    container->set_distributed_data(data);
	}
    }

    /*!
     * \brief Given a field, get a global data element to be be received from
     * a source.
     * \param field_name The name of the field to receive data from.
     * \param data The global data element.
     */
    void get_global_data(const std::string &field_name,
			 const DataType &data)
    {
	if ( field_name = "SCALAR_TEST_FIELD" )
	{
	    container->set_scalar_data(data);
	}
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// messenger test
void messenger_test(Parallel_Unit_Test &ut)
{
    // Create a data container instance.
    Teuchos::RCP<Data_Container> container = Teuchos::rcp(new Data_Container);

    // Create an instance of the source interface.
    Teuchos::RCP<Transfer_Data_Source<double> > tds = 
	Teuchos::rcp(new test_Transfer_Data_Source<double>());

    // Create an instance of the target interface.
    Teuchos::RCP<Transfer_Data_Target<double> > tdt = 
	Teuchos::rcp(new test_Transfer_Data_Target<double>(container));

    // Create a distributed field for these interfaces to be transferred.
    Teuchos::RCP<Transfer_Data_Field<double> > field = Teuchos::rcp(
	new Transfer_Data_Field<double>("DISTRIBUTED_TEST_FIELD", tds, tdt));

    // Create a map to populate.
    Teuchos::RCP<Transfer_Map> map = Teuchos::rcp(new Transfer_Map());

    // Create a mapper and populate the map.
    Mapper mapper;
    mapper.map(get_comm_world(), field, map);

    // Apply the map to the field.
    field->set_map(map);

    // Communicate the field with the messenger.
    Messenger messenger;
    messenger.communicate(get_comm_world(), field);

    // Check the results of the transfer.
    UNIT_TEST( container->get_distributed_handles().size() == 1);
    UNIT_TEST( container->get_distributed_handles[0] == nemesis::node() );
    UNIT_TEST( container->get_distributed_data().size() == 1);
    UNIT_TEST( container->get_distributed_data()[0] == 1.0*nemesis::node() );

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Messenger test ok on " << nemesis::node();
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, coupler::release);

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

	messenger_test(ut);
	gpass += ut.numPasses;
	gfail += ut.numFails;
	ut.reset();

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstMessenger, " 
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstMessenger, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstMessenger.cc
//---------------------------------------------------------------------------//
