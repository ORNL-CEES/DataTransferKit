//---------------------------------------------------------------------------//
/*!
 * \file DTK_DataTarget_def.hpp
 * \author Stuart R. Slattery
 * \brief DataTarget defintion.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATATARGET_DEF_HPP
#define DTK_DATATARGET_DEF_HPP

#include <DTK_Exception.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor.
 * \param mpi_comm Raw MPI communicator over which this data target is
 * defined. 
 */
template<class CoordinateField, class DataField>
DataTarget<CoordinateField,DataField>::DataTarget( const MPI_Comm& mpi_comm )
{
    // Wrap the raw communicator in a Teuchos comm object.
    Teuchos::RCP< Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	Teuchos::opaqueWrapper( mpi_comm );
    d_comm = Teuchos::rcp( new Teuchos::MpiComm<int>( opaque_comm ) );
    testPostcondition( d_comm != Teuchos::null,
		       "Error creating Teuchos comm from MPI comm." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class CoordinateField, class DataField>
DataTarget<CoordinateField,DataField>::~DataTarget()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Register a target field.
 * \param field_name The string key for this field.
 * \param coordinate_field The coordinate field over which data will be
 * evaluated. 
 * \param data_space The data space in which to write the evaluated data. This
 * space will neither be allocated or deallocated and should persist.
 */
template<class CoordinateField, class DataField>
void DataTarget<CoordinateField,DataField>::registerTargetField( 
    const std::string& field_name, 
    const CoordinateField& coordinate_field,
    DataField& data_space )
{
    // Create an integer id for this field.
    std::size_t field_id = d_name_map.size();
    
    // Add it to the maps.
    d_name_map[ field_name ] = field_id;
    d_coord_map[ field_id ] = coordinate_field;
    d_data_map[ field_id ] = data_space;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_DATATARGET_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_DataTarget_def.hpp
//---------------------------------------------------------------------------//
