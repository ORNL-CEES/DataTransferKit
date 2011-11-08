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

#include "Transfer_Map.hh"
#include "Message_Buffer.hh"
#include "Physics.hh"
#include "utils/SP.hh"
#include "utils/Packing_Utils.hh"
#include "comm/global.hh"

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
    //! Useful Typedefs.
    typedef DataType_T                          DataType;
    typedef int                                 OrdinateType;
    typedef int                                 HandleType;
    typedef Message_Buffer<OrdinateType>        Message_Buffer_t;
    typedef typename Message_Buffer_t::Buffer   Buffer;
    typedef denovo::SP<Physics>                 SP_Physics;
    typedef nemesis::Communicator_t             Communicator;
    //@}
    
  private:

    // Global communicator.
    const Communicator &d_comm_global;

    // Field name.
    const std::string &d_field_name;

    // Source physics.
    SP_Physics d_source;

    // Target physics.
    SP_Physics d_target;

 public:

    // Constructor.
    Messenger(const Communicator &comm_global,
	      const std::string &field_name,
	      SP_Physics source,
	      SP_Physics target);

    // Communicate the field from the source to the target.
    void communicate();

  private:

    // Private typedefs
    typedef std::list<Message_Buffer_t>    BufferList;

    // Post the receives
    void post_receives(BufferList &buffer_list);

    // Send the data
    void send(const KeyType& key);

    // Fill the map with the received data
    void fill_nodes(BufferList &buffer_list, const KeyType &key);
};

} // end namespace coupler

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#endif // core_Messenger_hh

//---------------------------------------------------------------------------//
//              end of core/Messenger.hh
//---------------------------------------------------------------------------//
