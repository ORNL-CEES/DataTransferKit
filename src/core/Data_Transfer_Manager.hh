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

#include "Transfer_Data_Field.hh"
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
    typedef Transfer_Data_Field<DataType>            Transfer_Data_Field_t;
    typedef denovo::SP<Transfer_Data_Field_t>        SP_Transfer_Data_Field;
    typedef denovo::SP<LG_Indexer>                   SP_LG_Indexer;
    typedef denovo::SP<Transfer_Map>                 SP_Transfer_Map;
    typedef nemesis::Communicator_t                  Communicator;
    //@}

  private:

    // Global communicator.
    const Communicator &d_comm_global;

 public:

    // Constructor.
    Data_Transfer_Manager(const Communicator &comm_global);

    // Destructor.
    ~Data_Transfer_Manager();

    // Transfer data associated with a distributed field from a source physics
    // to a target physics. 
    void distributed_transfer(SP_Transfer_Data_Field transfer_data_field);

    // Transfer a scalar field from a source physics to a target physics.
    void scalar_transfer(SP_Transfer_Data_Field transfer_data_field);
};

} // end namespace coupler

#endif // core_Data_Transfer_Manager_hh

//---------------------------------------------------------------------------//
//              end of core/Data_Transfer_Manager.hh
//---------------------------------------------------------------------------//
