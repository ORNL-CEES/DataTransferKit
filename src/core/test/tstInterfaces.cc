//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstInterfaces.cc
 * \author Stuart Slattery
 * \date   Thu Dec 01 16:50:04 2011
 * \brief  Unit tests for the data transfer pure virtual interfaces.
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
#include "../Transfer_Data_Source.hh"
#include "../Transfer_Data_Target.hh"

#include "Teuchos_RCP.hpp"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using coupler::Transfer_Data_Source;
using coupler::Transfer_Data_Target;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

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
    void register_comm(const Communicator &comm)
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

	if ( x > 0 && y > 0 && z > 0 )
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
	    std::vector<double> local_data(1, 1.0);
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
	    data = 1.0;
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

    double scalar_data;
    std::vector<double> received_data;
    std::vector<int> received_handles;

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
    void register_comm(const Communicator &comm)
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
	    std::vector<int> local_handles(1, 1);
	    std::vector<double> local_coords(3, 1.0);

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
	    received_handles = handles;
	    received_data = data;
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
	    scalar_data = data;
	}
    }

    // Test functions to determine whether the receive_data and get_global_data
    // methods acquired the correct data.
    std::vector<int> check_distributed_handles()
    {
	return received_handles;
    }

    std::vector<double> check_distributed_data()
    {
	return received_data;
    }

    double check_scalar_data()
    {
	return scalar_data;
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void source_interface_test(Parallel_Unit_Test &ut)
{
    // create an instance of the source interface.
    teuchos::RCP<Transfer_Data_Source<double> > source_iface = 
	new test_Transfer_Data_Source<double>();

    // test the interface methods.
    nemesis::Communicator_t source_comm;
    source_iface->register_comm(source_comm);
#ifdef MPI_COMM
    UNIT_TEST( source_comm == MPI_COMM_WORLD );
#else
    UNIT_TEST( source_comm == 1 );
#endif

    UNIT_TEST( source_iface->field_supported("DISTRIBUTED_TEST_FIELD") );
    UNIT_TEST( source_iface->field_supported("SCALAR_TEST_FIELD") );
    UNIT_TEST( !source_iface->field_supported("FOO_TEST_FIELD") );

    UNIT_TEST( source_iface->get_points(1, 1.0, 1.0, 1.0) );
    UNIT_TEST( !source_iface->get_points(1, -1.0, -1.0, -1.0) );

    std::vector<double> data_to_send;
    std::vector<int> handles_to_send(1, 1);
    source_iface->send_data("FOO_TEST_FIELD", handles_to_send, data_to_send);
    UNIT_TEST( data_to_send.size() == 0 );
    source_iface->send_data("DISTRIBUTED_TEST_FIELD", 
			    handles_to_send, data_to_send);
    UNIT_TEST( data_to_send.size() == 1);
    UNIT_TEST( data_to_send[0] = 1.0 );

    double global_scalar = 0.0;
    source_iface->set_global_data("FOO_TEST_FIELD", global_scalar);
    UNIT_TEST( global_scalar = 0.0 );
    source_iface->send_data("SCALAR_TEST_FIELD", global_scalar);
    UNIT_TEST( global_scalar = 1.0 );
    
    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Transfer_Data_Source test passes on " << node;
        ut.passes( m.str() );
    }
}

void target_interface_test(Parallel_Unit_Test &ut)
{
    // create an instance of the target interface.
    teuchos::RCP<Transfer_Data_Target<double> > target_iface = 
	new test_Transfer_Data_Target<double>();

    // test the interface methods.
    nemesis::Communicator_t target_comm;
    target_iface->register_comm(target_comm);
#ifdef COMM_MPI
    UNIT_TEST( target_comm == MPI_COMM_WORLD );
#else
    UNIT_TEST( target_comm == 1 );
#endif

    UNIT_TEST( target_iface->field_supported("DISTRIBUTED_TEST_FIELD") );
    UNIT_TEST( target_iface->field_supported("SCALAR_TEST_FIELD") );
    UNIT_TEST( !target_iface->field_supported("FOO_TEST_FIELD") );

    std::vector<double> target_coords;
    std::vector<int> target_handles;
    target_iface->set_points("FOO_TEST_FIELD", target_handles, target_coords);
    UNIT_TEST( target_coords.size() == 0 );
    UNIT_TEST( target_handles.size() == 0 );
    target_iface->set_points("DISTRIBUTED_TEST_FIELD", 
			     target_handles, target_coords);
    UNIT_TEST( target_coords.size() == 3 );
    UNIT_TEST( target_coords[0] == 1.0 );
    UNIT_TEST( target_coords[1] == 1.0 );
    UNIT_TEST( target_coords[2] == 1.0 );
    UNIT_TEST( target_handles.size() == 1 );
    UNIT_TEST( target_handles[0] == 1 );

    std::vector<double> data_to_receive(1, 1.0);
    std::vector<int> handles_to_receive(1, 1);
    target_iface->receive_data("FOO_TEST_FIELD", 
			       handles_to_receive, data_to_receive);
    UNIT_TEST( target_iface->check_distributed_data().size == 0 );
    UNIT_TEST( target_iface->check_distributed_handles().size == 0 );
    target_iface->receive_data("DISTRIBUTED_TEST_FIELD", 
			    handles_to_receive, data_to_receive);
    UNIT_TEST( target_iface->check_distributed_data().size == 1 );
    UNIT_TEST( target_iface->check_distributed_handles().size == 1 );
    UNIT_TEST( target_iface->check_distributed_data()[0] == 1.0 );
    UNIT_TEST( target_iface->check_distributed_handles()[0] == 1 );

    double global_scalar = 1.0;
    target_iface->get_global_data("FOO_TEST_FIELD", global_scalar);
    UNIT_TEST( target_iface->check_scalar_data() != 1.0 );
    target_iface->receive_data("SCALAR_TEST_FIELD", global_scalar);
    UNIT_TEST( target_iface->check_scalar_data() = 1.0 );

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Transfer_Data_Target test passes on " << node;
        ut.passes( m.str() );
    }
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, coupler::release);

    node  = nemesis::node();
    nodes = nemesis::nodes();
    
    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        source_interface_test(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

        target_interface_test(ut);
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
        std::cout << "ERROR: While testing tstInterfaces, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstInterfaces, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstInterfaces.cc
//---------------------------------------------------------------------------//
