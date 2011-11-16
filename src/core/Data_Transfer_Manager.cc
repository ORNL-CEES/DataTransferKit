//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Data_Transfer_Manager.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:44 2011
 * \brief  Data_Transfer_Manager member definitons.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Data_Transfer_Manager.hh"
#include "Mapper.hh"
#include "Messenger.hh"
#include "harness/DBC.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * \param comm_global Global communicator that encompasses all processes
 * participating in coupling. All methods driven by the manager will operate
 * on this communicator.
 */
template<class DataType_T>
Data_Transfer_Manager<DataType_T>::Data_Transfer_Manager(
    const Communicator &comm_global)
    : d_comm_global(comm_global)
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class DataType_T>
Data_Transfer_Manager<DataType_T>::~Data_Transfer_Manager()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Register a physics to be controlled by the manager.
 * \param physics_name Name of the physics being registered.
 * \param te Transfer_Evaluator implementation for the phyiscs being
 * registered.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::add_physics(
    const std::string &physics_name,
    SP_Transfer_Evaluator te)
{
    // Make a physics object.
    SP_Physics new_physics = new Physics<DataType>(te, d_comm_global);

    // Add it to the physics database.
    d_physics_db.insert( Physics_Pair(physics_name, new_physics) );
    
    Ensure( d_physics_db[physics_name] );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the topology map for transfer from a source physics to a
 * target physics for a particular field.
 * \param field_name the name of the field being mapped.
 * \param source_physics The name of the source physics used for the mapping.
 * \param target_physics The name of the target physics used for the mapping,
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::map(
    const std::string &field_name,
    const std::string &source_physics,
    const std::string &target_physics)
{
    // Get the physics that we are operating on.
    SP_Physics source = d_physics_db[source_physics];
    SP_Physics target = d_physics_db[target_physics];

    // Require that these physics support the field being mapped.
    Require( source->te()->field_supported(field_name) &&
	     target->te()->field_supported(field_name) );

    // Create a mapper.
    Mapper<DataType> mapper(d_comm_global, field_name, source, target);

    // Generate the map.
    mapper.map();

    Ensure( source->get_map(target_physics, field_name) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Transfer data associated with a field from a source physics to a
 * target physics.
 * \param field_name The name of the field being transferred.
 * \param source_physics The name of the source physics for the transfer.
 * \param target_physics The name of the target physics for the transfer.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::transfer(
    const std::string &field_name,
    const std::string &source_physics,
    const std::string &target_physics)
{
    // Get the physics we are operating on.
    SP_Physics source = d_physics_db[source_physics];
    SP_Physics target = d_physics_db[target_physics];

    // Require that these physics support the field being transferred.
    Require( source->te()->field_supported(field_name) &&
	     target->te()->field_supported(field_name) );
    
    // Create a messenger.
    Messenger<DataType> messenger(d_comm_global, field_name, source, target);

    // Transfer the field.
    messenger.communicate();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform a global rebalance on a field for conservation.
 * \param field_name The name of the field being balanced.
 * \param source_physics The name of the source physics.
 * \param target_physics The name of the target physics.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::balance(
    const std::string &field_name,
    const std::string &source_physics,
    const std::string &target_physics)
{
    // Get the physics we are operating on.
    SP_Physics source = d_physics_db[source_physics];
    SP_Physics target = d_physics_db[target_physics];

    // Require that these physics support the field being balanced.
    Require( source->te()->field_supported(field_name) &&
	     target->te()->field_supported(field_name) );

    // Get the normalization factor from the source.
    DataType norm;
    source->te()->integrate(field_name, norm);

    // Rebalance the target with the norm.
    target->te()->rebalance(field_name, norm);
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
