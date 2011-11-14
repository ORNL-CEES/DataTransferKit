//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Physics.cc
 * \author Stuart Slattery
 * \date   Wed Nov 02 11:01:32 2011
 * \brief  Physics member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Physics.hh"
#include "harness/DBC.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * \param physics_name The name of this physics.
 * \param te Transfer evaluator implementation for this physics.
 * \param comm_global The global communicator encapsulating all physics being
 * coupled. 
 */
Physics::Physics(const std::string &physics_name,
			     Transfer_Evaluator_t *te, 
			     Communicator comm_global)
    : d_name(physics_name)
{
    // Wrap the raw transfer evaluator pointer.
    d_te = te;

    // Get the sub communicator.
    d_te->register_comm(d_comm);

    // Generate local-to-global indexing.
    d_indexer = new LG_Indexer(comm_global, d_comm, d_te);
    Ensure(d_indexer);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Physics::~Physics()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Given a target physics and a field, add the mapping for which this
 * physics is the source. 
 * \param target_physics The name of the target physics being mapped to.
 * \param field_name The name of the field being mapped.
 * \param transfer_map Smart pointer to the Transfer_Map to be assigned to
 * this physics. 
 */
void Physics::set_map(std::string target_physics, 
		      std::string field_name,
		      SP_Transfer_Map transfer_map)
{
    if ( d_target_map[target_physics] )
    {
	d_target_map[target_physics].insert(
	    std::pair<std::string,SP_Transfer_Map>(field_name,transfer_map) );
    }

    else
    {
	Field_Map new_map;
	new_map.insert(
	    std::pair<std::string,SP_Transfer_Map>(field_name,transfer_map) );
	d_target_map.insert(
	    std::pair<std::string,Field_Map(target_physics,new_map) );
    }

    Ensure( (d_target_map[target_physics])[field_name] );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a target physics and a field, return the mapping for which
 * this physics is the source.
 * \param target_physics The name of the target physics.
 * \param field_name The name of the mapped field.
 * \return Returns a smart pointer to the Transfer_Map for the given target
 * physics and field.
 */
const Physics::SP_Transfer_Map Physics::get_map(std::string target_physics,
						std::string field_name)
{
    Require( (d_target_map[target_physics])[field_name] );
    return (d_target_map[target_physics])[field_name];
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Physics.cc
//---------------------------------------------------------------------------//
