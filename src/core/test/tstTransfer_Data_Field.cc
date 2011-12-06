//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstTransfer_Data_Field.cc
 * \author Stuart Slattery
 * \date   Fri Nov 18 14:43:10 2011
 * \brief  Transfer_Data_Field unit tests
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
#include "../Transfer_Map.hh"
#include "../Transfer_Data_Source.hh"
#include "../Transfer_Data_Target.hh"
#include "../Transfer_Data_Field.hh"

#include "Teuchos_RCP.hpp"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using coupler::Transfer_Map;
using coupler::Transfer_Data_Source;
using coupler::Transfer_Data_Target;
using coupler::Transfer_Data_Field;

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

    typedef double                                   DataType;
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;

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

    typedef double                                   DataType;
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;

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

    void set_points(const std::string &field_name,
		    std::vector<HandleType> &handles,
		    std::vector<CoordinateType> &coordinates)
    {
	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    std::vector<int> local_handles(1, 1);
	    std::vector<double> local_coords(3, 1.0);

	    handles = local_handles;
	    coordinates = local_coords;
	}
    }

    void receive_data(const std::string &field_name,
		      const std::vector<HandleType> &handles,
		      const std::vector<DataType> &data)
    {
	if ( field_name == "DISTRIBUTED_TEST_FIELD" )
	{
	    received_handles = handles;
	    received_data = data;
	}
    }

    void get_global_data(const std::string &field_name,
			 const DataType &data)
    {
	if ( field_name == "SCALAR_TEST_FIELD" )
	{
	    scalar_data = data;
	}
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void distributed_test(Parallel_Unit_Test &ut)
{
    // Create an instance of the source interface.
    Teuchos::RCP<Transfer_Data_Source<double> > tds = 
	Teuchos::rcp(new test_Transfer_Data_Source<double>());

    // Create an instance of the target interface.
    Teuchos::RCP<Transfer_Data_Target<double> > tdt = 
	Teuchos::rcp(new test_Transfer_Data_Target<double>());

    // Create a distributed field for these interfaces to be transferred.
    Transfer_Data_Field<double> field("DISTRIBUTED_TEST_FIELD", tds, tdt);

    // Add a transfer map to the field.
    UNIT_TEST( !field.is_mapped() );
    Teuchos::RCP<Transfer_Map> map = Teuchos::rcp(new Transfer_Map());
    field.set_map(map);

    // Test the functionality.
    UNIT_TEST( field.name() == "DISTRIBUTED_TEST_FIELD" );
    UNIT_TEST( field.source() == tds );
    UNIT_TEST( field.target() == tdt );
    UNIT_TEST( field.get_map() == map );
    UNIT_TEST( !field.is_scalar() );
    UNIT_TEST( field.is_mapped() );

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Distributed field test ok on " << nemesis::node();
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//
void scalar_test(Parallel_Unit_Test &ut)
{
    // Create an instance of the source interface.
    Teuchos::RCP<Transfer_Data_Source<double> > tds = 
	Teuchos::rcp(new test_Transfer_Data_Source<double>());

    // Create an instance of the target interface.
    Teuchos::RCP<Transfer_Data_Target<double> > tdt = 
	Teuchos::rcp(new test_Transfer_Data_Target<double>());

    // Create a scalar field for these interfaces to be transferred.
    Transfer_Data_Field<double> field("SCALAR_TEST_FIELD", tds, tdt, true);

    // Test the functionality.
    UNIT_TEST( field.name() == "SCALAR_TEST_FIELD" );
    UNIT_TEST( field.source() == tds );
    UNIT_TEST( field.target() == tdt );
    UNIT_TEST( field.is_scalar() );
    UNIT_TEST( !field.is_mapped() );

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Scalar field test ok on " << nemesis::node();
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
 
	distributed_test(ut);
	gpass += ut.numPasses;
	gfail += ut.numFails;
	ut.reset();

	scalar_test(ut);
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
        std::cout << "ERROR: While testing tstTransfer_Data_Field, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstTransfer_Data_Field, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstTransfer_Data_Field.cc
//---------------------------------------------------------------------------//
