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
// CONSTRUCTOR and DESTRUCTOR
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
// PUBLIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the topology map for transfer from a source physics to a
 * target physics for a particular field.
 * \param field_name the name of the field being mapped.
 * \param source The data transfer source implemenation.
 * \param target The data transfer target implemenation.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::map(
    const std::string &field_name,
    SP_Data_Transfer_Source source,
    SP_Data_Transfer_Target target)
{
    // Require that these physics support the field being mapped.
    Require( source->field_supported(field_name) &&
	     target->field_supported(field_name) );

    // Create a mapper.
    Mapper<DataType> mapper(d_comm_global, field_name, source, target);

    // Create a map.
    SP_Transfer_Map transfer_map = new Transfer_Map();

    // Generate the map.
    mapper.map(transfer_map);

    // Add the new map to the database.
    d_map_db->set_map(field_name, source, target, transfer_map);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Transfer data associated with a field from a source physics to a
 * target physics.
 * \param field_name The name of the field being transferred.
 * \param source The data transfer source implemenation.
 * \param target The data transfer target implemenation.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::distributed_transfer(
    const std::string &field_name,
    SP_Data_Transfer_Source source,
    SP_Data_Transfer_Target target)
{
    // Require that these physics support the field being transferred.
    Require( source->field_supported(field_name) &&
	     target->field_supported(field_name) );
    
    // Create a messenger.
    Messenger<DataType> messenger(d_comm_global, field_name, source, target);

    // Transfer the field.
    messenger.communicate( d_map_db->get_map(field_name, source, target) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform a global rebalance on a field for conservation.
 * \param field_name The name of the field being balanced.
 * \param source_physics The name of the source physics.
 * \param target_physics The name of the target physics.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::scalar_transfer(
    const std::string &field_name,
    SP_Data_Transfer_Source source,
    SP_Data_Transfer_Target target)
{
    // Require that these physics support the field being balanced.
    Require( source->field_supported(field_name) &&
	     target->field_supported(field_name) );

    // Source sets the scalar.
    DataType data;
    source->set_global_data(field_name, data);

    // Target gets the scalar.
    target->get_global_data(field_name, data);
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Given a target physics and a field, add the mapping for which this
 * physics is the source. 
 * \param field_name The name of the mapped field.
 * \param source The data transfer source implemenation.
 * \param target The data transfer target implemenation.
 * \param transfer_map Smart pointer to the Transfer_Map created for the
 * transfer. 
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::set_map(const std::string &field_name,
						SP_Data_Transfer_Source source,
						SP_Data_Transfer_Target target,
						SP_Transfer_Map transfer_map)
{
    if ( d_map_db[target_physics] )
    {
	d_map_db[target_physics].insert(
	    std::pair<std::string,SP_Transfer_Map>(field_name,transfer_map) );
    }

    else
    {
	Field_Map new_map;
	new_map.insert(
	    std::pair<std::string,SP_Transfer_Map>(field_name,transfer_map) );
	d_map_db.insert(
	    std::pair<std::string,Field_Map>(target_physics,new_map) );
    }

    Ensure( (d_map_db[target_physics])[field_name] );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a source, target and field, return the mapping.
 * \param field_name The name of the mapped field.
 * \param source The data transfer source implemenation.
 * \param target The data transfer target implemenation.
 * \return Returns a smart pointer to a Transfer_Map.
 */
template<class DataType_T>
const SP_Transfer_Map 
Data_Transfer_Manager<DataType_T>::get_map(const std::string &field_name,
					   SP_Data_Transfer_Source source,
					   SP_Data_Transfer_Target target)
{
    Require( d_map_db[field_name][source][target] );
    return d_map_db[field_name][source][target];
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
