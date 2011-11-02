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

#include <vector>
#include <map>
#include <string>

#include "Transfer_Evaluator.hh"
#include "Transfer_Map.hh"
#include "Physics.hh"
#include "LG_Indexer.hh"
#include "Messenger.hh"
#include "field/Field.hh"
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
    typedef FieldType_T                              FieldType;
    typedef int                                      Handle;
    typedef const Handle*                            Handle_Iterator;
    typedef double                                   Coordinate;
    typedef const Coordinate*                        Coord_Iterator;
    typedef denovo::SP<Physics>                      SP_Physics;
    typedef denovo::SP<Messenger<FieldType> >        SP_Messenger;
    typedef denovo::SP<Transfer_Map>                 SP_Transfer_Map;
    typedef nemesis::Communicator_t                  Communicator;
    typedef std::vector<char>                        Buffer;
    typedef std::map<std::string,SP_Physics>         Physics_DB;
    //@}

  private:

    // Global communicator.
    Communicator d_comm_global;

    // Physics database.
    Physics_DB d_physics_db;
    
  public:

    //! Constructor.
    Data_Transfer_Manager(Communicator comm_global);

    //! Destructor.
    ~Data_Transfer_Manager();

    //! Register a physics with the manager.
    void add_physics(std::string physics_name, 
		     Transfer_Evaluator *te);

    //! Register a field with the manager.
    void add_field(std::string field_name);

    //! Build the topology map for transfer from a source physics to a target
    //! physics for a particular field.
    void map(std::string field_name,
	     std::string source_physics,
	     std::string target_physics);

    //! Transfer data associated with a field from a source physics to a target
    //! physics. 
    void transfer(std::string field_name,
		  std::string source_physics,
		  std::string target_physics);
};

} // end namespace coupler

#endif // coupler_Data_Transfer_Manager_hh

//---------------------------------------------------------------------------//
//              end of core/Data_Transfer_Manager.hh
//---------------------------------------------------------------------------//
