//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Data_Transfer_Manager.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:44 2011
 * \brief  Data_Transfer_Manager class definiton.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Data_Transfer_Manager_hh
#define core_Data_Transfer_Manager_hh

#include <vector>
#include <map>
#include <string>

#include "Transfer_Data_Source.hh"
#include "Transfer_Data_Target.hh"
#include "Transfer_Map.hh"
#include "LG_Indexer.hh"
#include "utils/SP.hh"
#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Data_Transfer_Manager
 * \brief The Data_Transfer_Manager manages the data transfer problem between
 * physics. 
 *
 * All objects required for coupling are organized by the manager. The manager
 * is the top level mechanism for interacting with the coupling package.
 */
//===========================================================================//

template<class DataType_T>
class Data_Transfer_Manager 
{
  public:
    
    //@{
    //! Useful typedefs.
    typedef DataType_T                               DataType;
    typedef Data_Transfer_Source<DataType>           Data_Transfer_Source_t;
    typedef denovo::SP<Data_Transfer_Source_t>       SP_Data_Transfer_Source;
    typedef Data_Transfer_Target<DataType>           Data_Transfer_Target_t;
    typedef denovo::SP<Data_Transfer_Target_t>       SP_Data_Transfer_Target;
    typedef denovo::SP<LG_Indexer>                   SP_LG_Indexer;
    typedef denovo::SP<Transfer_Map>                 SP_Transfer_Map;
    typedef std::map<SP_Data_Transfer_Target,SP_Transfer_Map>   Target_Map;
    typedef std::map<SP_Data_Transfer_Source,Target_Map>        Source_Map;
    typedef std::map<std::string,Source_Map>         Map_DB;
    typedef nemesis::Communicator_t                  Communicator;
    //@}

  private:

    // Global communicator.
    const Communicator &d_comm_global;

    // Map database.
    Map_DB d_map_db;
    
  public:

    // Constructor.
    Data_Transfer_Manager(const Communicator &comm_global);

    // Destructor.
    ~Data_Transfer_Manager();

    // Build the topology map for transfer from a source physics to a target
    // physics for a particular field.
    void map(const std::string &field_name,
	     SP_Data_Transfer_Source source,
	     SP_Data_Transfer_Target target);

    // Transfer data associated with a distributed field from a source physics
    // to a target physics. 
    void distributed_transfer(const std::string &field_name,
			      SP_Data_Transfer_Source source,
			      SP_Data_Transfer_Target target);

    // Transfer a scalar field from a source physics to a target physics.
    void scalar_transfer(const std::string &field_name,
			 SP_Data_Transfer_Source source,
			 SP_Data_Transfer_Target target); 

  private:

    // Given a target physics and a field, add the mapping for which this
    // physics is the source.
    void set_map(const std::string &field_name,
		 SP_Data_Transfer_Source source,
		 SP_Data_Transfer_Target target,
		 SP_Transfer_Map transfer_map);

    // Given a target physics and a field, return the mapping for which this
    // physics is the source.
    const SP_Transfer_Map get_map(const std::string &field_name,
				  SP_Data_Transfer_Source source,
				  SP_Data_Transfer_Target target);
};

} // end namespace coupler

#endif // core_Data_Transfer_Manager_hh

//---------------------------------------------------------------------------//
//              end of core/Data_Transfer_Manager.hh
//---------------------------------------------------------------------------//
