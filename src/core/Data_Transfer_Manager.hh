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

#include "Transfer_Evaluator.hh"
#include "Physics.hh"
#include "Mapper.hh"
#include "Messenger.hh"
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
    typedef Transfer_Evaluator<DataType>             Transfer_Evaluator_t;
    typedef Physics<DataType>                        Physics_t;
    typedef denovo::SP<Physics_t>                    SP_Physics;
    typedef std::pair<std::string,SP_Physics>        Physics_Pair;
    typedef std::map<std::string,SP_Physics>         Physics_DB;
    typedef nemesis::Communicator_t                  Communicator;
    //@}

  private:

    // Global communicator.
    const Communicator &d_comm_global;

    // Physics database.
    Physics_DB d_physics_db;
    
  public:

    //! Constructor.
    Data_Transfer_Manager(const Communicator &comm_global);

    //! Destructor.
    ~Data_Transfer_Manager();

    //! Register a physics with the manager.
    void add_physics(const std::string &physics_name, 
		     Transfer_Evaluator_t *te);

    //! Build the topology map for transfer from a source physics to a target
    //! physics for a particular field.
    void map(const std::string &field_name,
	     const std::string &source_physics,
	     const std::string &target_physics);

    //! Transfer data associated with a field from a source physics to a
    //! target physics. 
    void transfer(const std::string &field_name,
		  const std::string &source_physics,
		  const std::string &target_physics);
};

} // end namespace coupler

#endif // core_Data_Transfer_Manager_hh

//---------------------------------------------------------------------------//
//              end of core/Data_Transfer_Manager.hh
//---------------------------------------------------------------------------//
