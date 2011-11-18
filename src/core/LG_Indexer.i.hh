//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/LG_Indexer.i.hh
 * \author Stuart R. Slattery
 * \date   Thu Jun 16 16:23:46 2011
 * \brief  Member definitions of class LG_Indexer.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.i.hh,v 1.4 2008/01/04 22:50:12 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_LG_Indexer_i_hh
#define core_LG_Indexer_i_hh

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief LG_Indexer constructor.
 *
 * \param comm_global Communicator that contains the entire coupled simulation.
 * \param comm_local Communicator for the local application being indexed.
 */
template<class LocalApp>
LG_Indexer::LG_Indexer(const Communicator &comm_global, 
                       const Communicator &comm_local,
                       denovo::SP<LocalApp> local_app)
{
    // Indicate whether we have the local app.
    int local_app_indicator = 0;
    if ( local_app )
    {
        local_app_indicator = 1;
    }

    // Get the local id from comm_local.
    nemesis::set_internal_comm(comm_local);
    OrdinateType local_id = nemesis::node();

    // Get global information from comm_global.
    nemesis::set_internal_comm(comm_global);

    // Make a vector of local ids everywhere.
    Vec_Ordinate local_ids(nemesis::nodes(), 0);
    local_ids[nemesis::node()] = local_id;
    nemesis::global_sum(&local_ids[0], nemesis::nodes());

    // Make a vector of application indicators.
    Vec_Ordinate app_ids(nemesis::nodes(), 0);
    app_ids[nemesis::node()] = local_app_indicator;
    nemesis::global_sum(&app_ids[0], nemesis::nodes());

    // Make a vector of global ids everywhere.
    Vec_Ordinate global_ids(nemesis::nodes(), 0);
    global_ids[nemesis::node()] = nemesis::node();
    nemesis::global_sum(&global_ids[0], nemesis::nodes());

    Check( local_ids.size() == global_ids.size() );
    Check( app_ids.size() == global_ids.size() );

    // Make the map.
    Ordinate_Iterator app_iter = app_ids.begin();
    Ordinate_Iterator global_iter = global_ids.begin();
    for (Ordinate_Iterator local_iter = local_ids.begin(), 
		       local_iter_end = local_ids.end();
         local_iter != local_iter_end; 
	 ++local_iter)
    {
        if (*app_iter)
        {
            d_l2g_map[*local_iter] = *global_iter;
        }

        ++global_iter, ++app_iter;
    }
}

} // end namespace coupler

#endif // core_LG_Indexer_i_hh

//---------------------------------------------------------------------------//
//              end of core/LG_Indexer.i.hh
//---------------------------------------------------------------------------//
