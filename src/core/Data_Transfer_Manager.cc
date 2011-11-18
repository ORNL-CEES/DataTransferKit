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
 * \brief Transfer data associated with a field from a source physics to a
 * target physics.
 * \param transfer_data_field The field being transfered.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::distributed_transfer(
    SP_Transfer_Data_Field transfer_data_field)
{
    // Require that the field is destributed.
    Require( !transfer_data_field->is_scalar() );

    // If the field has not yet been mapped from the source to the target,
    // perform the mapping before the transfer.
    if ( !transfer_data_field->is_mapped() )
    {
	// Create a mapper.
	Mapper<DataType> mapper;

	// Create a map.
	SP_Transfer_Map transfer_map = new Transfer_Map();

	// Generate the map.
	mapper.map(d_comm_global, transfer_data_field, transfer_map);

	// Add the new map to the database.
	transfer_data_field->set_map(transfer_map);
    }
    
    // Check that the mapping has been applied.
    Check ( transfer_data_field->is_mapped() )

    // Create a messenger.
    Messenger<DataType> messenger(d_comm_global, transfer_data_field);

    // Transfer the field.
    messenger.communicate();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send a global scalar quantity from the source to the target.
 * \param transfer_data_field The the scalar field being transfered.
 */
template<class DataType_T>
void Data_Transfer_Manager<DataType_T>::scalar_transfer(
    SP_Transfer_Data_Field transfer_data_field)
{
    // Require that the field is scalar.
    Require( transfer_data_field->is_scalar() );

    // Source sets the scalar.
    DataType data;
    transfer_data_field->source()->set_global_data(field_name, data);

    // Target gets the scalar.
    transfer_data_field->target()->get_global_data(field_name, data);
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
