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
#include "Physics.hh"
#include "LG_Indexer.hh"
#include "Messenger.hh"
#include "field/Field.hh"
#include "fields/Field_DB.hh"
#include "utils/SP.hh"
#include "comm/global.hh"
#include <vector>
#include <map>
#include <string>

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
    typedef FieldType_T                              FieldType;
    typedef int                                      Handle;
    typedef const Handle*                            Handle_Iterator;
    typedef double                                   Coordinate;
    typedef const Coordinate*                        Coord_Iterator;
    typedef denovo::SP<Transfer_Evaluator>           SP_Transfer_Evaluator;
    typedef denovo::SP<Physics>                      SP_Physics;
    typedef Field_DB<FieldType>                      DB;
    typedef denovo::SP<DB>                           SP_DB;
    typedef denovo::SP<LG_Indexer>                   SP_LG_Indexer;
    typedef denovo::SP<Messenger>                    SP_Messenger;
    typedef nemesis::Communicator_t                  Communicator;
    typedef std::vector<char>                        Buffer;
    typedef std::map<std::string,SP_Physics>         Physics_DB;
    //@}

  private:

    // Global communicator.
    Communicator d_comm_global;

    // Physics A communicator.
    Communicator d_comm_a;
   
    // Physics B communicator.
    Communicator d_comm_b;

    // Physics A transfer evaluator.
    SP_Transfer_Evaluator d_te_a;

    // Physics B transfer evaluator.
    SP_Transfer_Evaluator d_te_b;

    // Physics A local to global indexer.
    SP_LG_Indexer d_indexer_a;

    // Physics B local to global indexer.
    SP_LG_Indexer d_indexer_b;

    // Physics A messenger object.
    SP_Messenger d_messenger_a;

    // Physics B messenger object.
    SP_Messenger d_messenger_b;

    // Field database.
    SP_DB d_f_db;

    // Physics database.
    Physics_DB d_physics_db;
    

  public:

    // Constructor.
    Data_Transfer_Manager(Communicator comm_global);

    // Destructor.
    ~Data_Transfer_Manager();

    // Register a physics with the manager.
    void add_physics(std::string physics_name, 
		     Transfer_Evaluator* te);

    // Register a field with the manager.
    void add_field(std::string field_name);

    // Build the topology map for transfer from a source physics to a target
    // physics.
    void map(std::string field_name,
	     std::string source_physics,
	     std::string target_physics);

    // Transfer data from a source physics to a target physics.
    void transfer(std::string field_name,
		  std::string source_physics,
		  std::string target_physics);
};

} // end namespace coupler

#endif // coupler_Data_Transfer_Manager_hh

//---------------------------------------------------------------------------//
//              end of core/Data_Transfer_Manager.hh
//---------------------------------------------------------------------------//
