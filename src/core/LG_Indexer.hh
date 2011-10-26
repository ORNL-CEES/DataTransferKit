//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/LG_Indexer.hh
 * \author Stuart R. Slattery
 * \date   Thu Jun 16 16:23:46 2011
 * \brief  LG_Indexer class definition
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_LG_Indexer_hh
#define coupler_LG_Indexer_hh

#include <map>
#include "utils/SP.hh"
#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class LG_Indexer
 * \brief Convert local-to-application process ID's to global process ID's
 *
 */
/*! 
 * \example core/test/tstLG_Indexer.cc
 *
 */
//===========================================================================//

class LG_Indexer 
{
  public:
    //@{
    //! Useful typedefs.
    typedef nemesis::Communicator_t     Communicator_t;
    //@}

  private:
    // Private typedef.
    typedef std::map<int, int>          Indexer_Map;

    // Local to global map
    Indexer_Map d_l2g_map;

  public:
    //! Constructor.
    template<class LocalApp>
    LG_Indexer(const Communicator_t &comm_world, 
               const Communicator_t &comm_local,
               denovo::SP<LocalApp> local_app);

    //! Convert application-local PID to global PID.
    int l2g(int local_pid) { return d_l2g_map[local_pid]; }

    //! Get the size of the indexer.
    int size() { return d_l2g_map.size(); }
};

} // end namespace coupler

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS AND TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "LG_Indexer.i.hh"

#endif // coupler_LG_Indexer_hh

//---------------------------------------------------------------------------//
//              end of coupler/LG_Indexer.hh
//---------------------------------------------------------------------------//
