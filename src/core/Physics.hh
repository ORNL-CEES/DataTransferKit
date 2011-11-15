//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Physics.hh
 * \author Stuart Slattery
 * \date   Wed Nov 02 11:01:32 2011
 * \brief  Physics class definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Physics_hh
#define core_Physics_hh

#include <string>
#include <map>

#include "Transfer_Evaluator.hh"
#include "Transfer_Map.hh"
#include "LG_Indexer.hh"
#include "utils/SP.hh"
#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Physics
 * \brief The Physics class is a container for all of the objects created for
 * each physics that participates in coupling. In addition, it serves as a
 * database for the mappings generated with the particular physics.
 *
 */
//===========================================================================//

template<class DataType_T>
class Physics 
{
  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                               DataType;
    typedef nemesis::Communicator_t                  Communicator;
    typedef Transfer_Evaluator<DataType>             Transfer_Evaluator_t;
    typedef denovo::SP<Transfer_Evaluator_t>         SP_Transfer_Evaluator;
    typedef denovo::SP<LG_Indexer>                   SP_LG_Indexer;
    typedef denovo::SP<Transfer_Map>                 SP_Transfer_Map;
    typedef std::map<std::string,SP_Transfer_Map>    Field_Map;
    typedef std::map<std::string,Field_Map>          Target_Map;
    //@}

  private:

    // Physics name.
    std::string d_name;

    // Transfer evaluator implementation.
    SP_Transfer_Evaluator d_te;

    // Sub-communicator.
    Communicator d_comm;

    // Local-to-global indexer.
    SP_LG_Indexer d_indexer;

    // Map database. This physics is the source. The key value designates the
    // target physics. The return value provides another database which is in
    // turn queried by using the field as the key value. This returns a
    // transfer map.
    Target_Map d_target_map;

  public:

    // Constructor.
    Physics(std::string physics_name, 
	    SP_Transfer_Evaluator te, 
	    Communicator comm_global);

    // Destructor.
    ~Physics();

    //! Return the name of this physics.
    const std::string& name() { return d_name; }

    //! Return the transfer evaluator implementation.
    const SP_Transfer_Evaluator te() { return d_te; }

    //! Return the communicator.
    const Communicator& comm() { return d_comm; }

    //! Return the indexer.
    const SP_LG_Indexer& indexer() { return d_indexer; }

    // Given a target physics and a field, add the mapping for which this
    // physics is the source.
    void set_map(std::string target_physics, 
		 std::string field_name,
		 SP_Transfer_Map transfer_map);

    // Given a target physics and a field, return the mapping for which this
    // physics is the source.
    const SP_Transfer_Map get_map(std::string target_physics,
				  std::string field_name);
};

} // end namespace coupler

#endif // core_Physics_hh

//---------------------------------------------------------------------------//
//              end of core/Physics.hh
//---------------------------------------------------------------------------//
