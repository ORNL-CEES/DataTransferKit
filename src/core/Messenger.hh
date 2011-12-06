//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Messenger.hh
 * \author Stuart R. Slattery
 * \date   Thu May 26 11:02:57 2011
 * \brief  Messenger class definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Messenger_hh
#define core_Messenger_hh

#include <string>
#include <list>

#include "Transfer_Map.hh"
#include "Transfer_Data_Field.hh"
#include "Message_Buffer.hh"
#include "Packing_Utils.hh"
#include "comm/global.hh"

#include "Teuchos_RCP.hpp"

namespace coupler
{

//===========================================================================//
/*!
 * \class Messenger
 * \brief A Messenger class to send and receive field data in a map
 *
 * \sa Messenger.cc for detailed descriptions.
 */
/*! 
 * \example coupler/test/tstMessenger.cc
 *
 * Test of Messenger.
 */
//===========================================================================//

template<class DataType_T>
class Messenger 
{
  public:   

    //@{
    //! Useful typedefs.
    typedef DataType_T                               DataType;
    typedef int                                      HandleType;
    typedef int                                      OrdinateType;
    typedef Message_Buffer<OrdinateType>             Message_Buffer_t;
    typedef typename Message_Buffer_t::Buffer        Buffer;
    typedef std::list<Message_Buffer_t>              BufferList;
    typedef typename BufferList::iterator            BufferList_Iterator;
    typedef Transfer_Data_Field<DataType>            Transfer_Data_Field_t;
    typedef Teuchos::RCP<Transfer_Data_Field_t>      RCP_Transfer_Data_Field;
    typedef Teuchos::RCP<Transfer_Map>               RCP_Transfer_Map;
    typedef typename Transfer_Map::Map_Iterator      Map_Iterator;
    typedef typename Transfer_Map::Map_Pair          Map_Pair;
    typedef typename Transfer_Map::Set_Iterator      Set_Iterator;
    typedef typename Transfer_Map::Set_Pair          Set_Pair;
    typedef nemesis::Communicator_t                  Communicator;
    //@}
    
 public:

    // Constructor.
    Messenger();

    // Destructor.
    ~Messenger();

    // Communicate the field from the source to the target.
    void communicate(const Communicator &comm_global,
		     RCP_Transfer_Data_Field transfer_data_field);

  private:

    // Target post receives.
    void post_receives(RCP_Transfer_Data_Field transfer_data_field,
		       BufferList &buffer_list);

    // Source send data.
    void send(RCP_Transfer_Data_Field transfer_data_field);

    // Process the target requests.
    void process_requests(RCP_Transfer_Data_Field transfer_data_field,
			  BufferList &buffer_list);
};

} // end namespace coupler

#endif // core_Messenger_hh

//---------------------------------------------------------------------------//
//              end of core/Messenger.hh
//---------------------------------------------------------------------------//
