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

#ifndef coupler_Data_Transfer_Manager_hh
#define coupler_Data_Transfer_Manager_hh

#include "Transfer_Evaluator.hh"
#include "Transfer_Map.hh"
#include "LG_Indexer.hh"
#include "Messenger.hh"
#include "fields/Field_DB.hh"
#include "utils/SP.hh"
#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Data_Transfer_Manager
 * \brief The Data_Transfer_Manager manages the data transfer problem between
 * codes. 
 *
 * All objects required for coupling are organized by the manager. The manager
 * is the top level mechanism for interacting with the coupling package.
 */
//===========================================================================//

template<class FieldType_T>
class Data_Transfer_Manager 
{
  public:
    
    //@{
    //! Useful typedefs.
    typedef Field_Type_T                             FieldType;
    typedef int                                      Handle;
    typedef const Handle*                            Handle_Iterator;
    typedef double                                   Coordinate;
    typedef const Coordinate*                        Coord_Iterator;
    typedef nemesis::SP<Transfer_Evaluator>          SP_Transfer_Evaluator;
    typedef Transfer_Map<Handle,int>                 Map;
    typedef nemesis::SP<Map>                         SP_Transfer_Map;
    typedef Field_DB<FieldType>                      DB;
    typedef nemesis::SP<DB>                          SP_DB;
    typedef nemesis::SP<LG_Indexer>                  SP_LG_Indexer;
    typedef nemesis::SP<Messenger>                   SP_Messenger;
    typedef nemesis::Communicator_t                  Communicator;
    typedef std::vector<char>                        Buffer;
    //@}

  private:

    // Global communicator.
    Communicator_t d_comm_global;

    // Physics A communicator.
    Communicator_t d_comm_a;
   
    // Physics B communicator.
    Communicator_t d_comm_b;

    // Physics A transfer evaluator.
    SP_Transfer_Evaluator d_te_a;

    // Physics B transfer evaluator.
    SP_Transfer_Evaluator d_te_b;

    // Physics A local to global indexer.
    SP_LG_Indexer d_indexer_a;

    // Physics B local to global indexer.
    SP_LG_Indexer d_indexer_b;

    // Topology map for transfer from A to B.
    SP_Transfer_Map d_map_a2b;

    // Topology map for transfer from B to A.
    SP_Transfer_Map d_map_b2a;

    // Physics A messenger object.
    SP_Messenger d_messenger_a;

    // Physics B messenger object.
    SP_Messenger d_messenger_b;

    // Field database.
    SP_DB d_f_db;

  public:

    // Constructor.
    Data_Transfer_Manager(Communicator_t comm_global,
	                  Transfer_Evaluator* TE_A,
			  Transfer_Evaluator* TE_B);

    // Destructor.
    ~Data_Transfer_Manager();

    // Register a field to the manager.
    void add_field(std::string field_name);

    // Build the topology map for transfer from A to B.
    void map_A2B(std::string field_name);

    // Build the topology map for transfer from B to A.
    void map_B2A(std::string field_name);

    // Transfer data from A to B.
    void transfer_A2B(std::string field_name);

    // Transfer data from B to A.
    void transfer_B2A(std::string field_name);
};

} // end namespace coupler

#endif // coupler_Data_Transfer_Manager_hh

//---------------------------------------------------------------------------//
//              end of core/Data_Transfer_Manager.hh
//---------------------------------------------------------------------------//
