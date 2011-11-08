//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Mapper.hh
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:31:05 2011
 * \brief  Mapper class defintion.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Mapper_hh
#define core_Mapper_hh

#include <string>
#include <list>

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
 * \class Mapper
 * \brief A mapper class to generate a one way mapping for transfer between a
 * source physics and a target physics.
 *
 * \sa Mapper.cc for detailed descriptions.
 */
/*! 
 * \example core/test/tstMapper.cc
 *
 * Test of Mapper.
 */
//===========================================================================//

class Mapper 
{
  public:

    //@{
    //! Useful typedefs.
    typedef int                                   HandleType;
    typedef const HandleType*                     Handle_Iterator;
    typedef int                                   OrdinateType;
    typedef double                                CoordinateType;
    typedef const CoordinateType*                 Coord_Iterator;
    typedef Message_Buffer<OrdinateType>          Message_Buffer_t;
    typedef typename Message_Buffer_t::Buffer     Buffer;
    typedef std::list<Message_Buffer_t>           BufferList;
    typedef typename BufferList::iterator         BufferList_Iterator;
    typedef denovo::SP<Physics>                   SP_Physics;
    typedef nemesis::Communicator_t               Communicator;
    typedef typename Transfer_Map::Map_Iterator   Map_Iterator;
    typedef typename Transfer_Map::Map_Pair       Map_Pair;
    typedef typename Transfer_Map::Set_Iterator   Set_Iterator;
    typedef typename Transfer_Map::Set_Pair       Set_Pair;
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
    Mapper(const Communicator &comm_global,
	   const std::string &field_name,
	   SP_Physics source,
	   SP_Physics target);

    // Map the field from the source onto the physics.
    void map();

  private:

    

};

} // end namespace coupler

#endif // core_Mapper_hh

//---------------------------------------------------------------------------//
//              end of core/Mapper.hh
//---------------------------------------------------------------------------//
