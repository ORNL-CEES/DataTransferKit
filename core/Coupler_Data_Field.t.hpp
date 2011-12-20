//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Coupler_Data_Field.t.hpp
 * \author Stuart Slattery
 * \date   Fri Nov 18 11:57:58 2011
 * \brief  Coupler_Data_Field template member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef core_Coupler_Data_Field_t_hpp
#define core_Coupler_Data_Field_t_hpp

#include <cassert>
#include <algorithm>
#include <vector>

#include <Mesh_SerializationTraits.hpp>

#include "Teuchos_CommHelpers.hpp"

#include "Tpetra_Vector.hpp"

namespace Coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR and DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * \param field_name The name of this field. Required by the
 * Coupler_Data_Source and Coupler_Data_Target interfaces to check field
 * support.
 * \param source Coupler_Data_Source implementation that will serve as the
 * data source for this field.
 * \param target Coupler_Data_Target implementation that will serve as the
 * target for this field.
 * \param scalar Set to true if this field is scalar, false if distributed.
 */
template<class DataType, class HandleType, class CoordinateType>
Data_Field<DataType,HandleType,CoordinateType>::Data_Field(
    RCP_Communicator comm_global,
    const std::string &field_name,
    RCP_Data_Source source,
    RCP_Data_Target target,
    bool scalar)
    : d_comm(comm_global)
    , d_field_name(field_name)
    , d_source(source)
    , d_target(target)
    , d_scalar(scalar)
    , d_mapped(false)
{ 
    // Require that these physics support the field being mapped.
    assert( d_source->field_supported(d_field_name) &&
	    d_target->field_supported(d_field_name) );

    // If not a scalar field, create the mapping between the data source and
    // data target.
    if ( !d_scalar )
    {
	point_map();
	d_mapped = true;
    }
}

/*!
 * \brief Destructor.
 */
template<class DataType, class HandleType, class CoordinateType>
Data_Field<DataType,HandleType,CoordinateType>::~Data_Field()
{ /* ... */ }

//---------------------------------------------------------------------------//
// PUBLIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Transfer data from the data source to the data target.
 */
template<class DataType, class HandleType, class CoordinateType>
void Data_Field<DataType,HandleType,CoordinateType>::transfer()
{ 
    // Scalar transfer.
    if ( d_scalar )
    {
	scalar_transfer();
    }

    // Distributed transfer.
    else
    {
	assert( d_mapped );
	distributed_transfer();
    }
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Generate topology map for this field based on point mapping.
 */
template<class DataType, class HandleType, class CoordinateType>
void Data_Field<DataType,HandleType,CoordinateType>::point_map()
{
    // Extract the local list of handles. These are the global indices for the
    // Tpetra map.
    const Teuchos::ArrayView<PointType> target_points = 
	d_target->set_points( d_field_name );
    typename Teuchos::ArrayView<PointType>::const_iterator target_point_it;

    std::vector<HandleType> target_handles(target_points.size());
    typename std::vector<HandleType>::iterator target_handle_it;

    for (target_handle_it = target_handles.begin(), 
	  target_point_it = target_points.begin(); 
	 target_handle_it != target_handles.end();
	 ++target_handle_it, ++target_point_it)
    {
	*target_handle_it = target_point_it->handle();
    }

    // Generate the map for the data target.
    const Teuchos::ArrayView<const HandleType> 	target_handles_view(target_handles);
    d_target_map = 
	Tpetra::createNonContigMap<HandleType>( target_handles_view, d_comm);

    // Data target communicate points to the data source.
    int local_size = target_points.size();
    int global_max = 0;
    Teuchos::reduceAll<OrdinalType,int>(*d_comm,
					Teuchos::REDUCE_MAX, 
					int(1), 
					&local_size, 
					&global_max);

    // Generate a target point buffer to send to the source. Pad the rest of
    // the buffer with null points. This is using -1 as the handle for a
    // null point here. This is OK as tpetra expects ordinals to be equal to
    // or greater than 0.
    PointType null_point(-1, 0.0, 0.0, 0.0);
    std::vector<PointType> send_points(global_max, null_point);
    typename std::vector<PointType>::iterator send_point_it;
    for (send_point_it = send_points.begin(),
	target_point_it = target_points.begin();
	 target_point_it != target_points.end();
	 ++send_point_it, ++target_point_it)
    {
	*send_point_it = *target_point_it;
    }

    // Communicate local points to all processes.
    std::vector<HandleType> source_handles;
    std::vector<PointType> receive_points(global_max, null_point);
    typename std::vector<PointType>::const_iterator receive_point_it;
    for ( int i = 0; i < d_comm->getSize(); ++i )
    {
	// Process rank i is broadcasting.
	if ( d_comm->getRank() == i )
	{
	    receive_points = send_points;
	}

	// Barrier before moving on to the next broadcast rank.
	Teuchos::barrier<OrdinalType>(*d_comm);

	Teuchos::broadcast<OrdinalType,PointType>( *d_comm,
						   i,
						   global_max,
						   &receive_points[0]);
						   
	// The data source searches for the points in its domain.
	for ( receive_point_it = receive_points.begin();
	      receive_point_it != receive_points.end();
	      ++receive_point_it )
	{
	    if ( receive_point_it->handle() > -1 )
	    {
		if ( d_source->get_points(*receive_point_it) )
		{
		    source_handles.push_back( receive_point_it->handle() );
		}
	    }
	}
    }

    // Barrier again before map generation.
    Teuchos::barrier<OrdinalType>(*d_comm);

    // Generate the map for the data source.
    const Teuchos::ArrayView<const HandleType> 
	source_handles_view(source_handles);

    d_source_map = 
	Tpetra::createNonContigMap<HandleType>( source_handles_view, d_comm);

    // Setup the exporter.
    d_export = Teuchos::rcp( 
	new Tpetra::Export<HandleType>(d_source_map, d_target_map) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform scalar transfer.
 */
template<class DataType, class HandleType, class CoordinateType>
void Data_Field<DataType,HandleType,CoordinateType>::scalar_transfer()
{
    d_target->get_global_data( d_field_name, 
			       d_source->set_global_data( d_field_name) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform distributed transfer.
 */
template<class DataType, class HandleType, class CoordinateType>
void Data_Field<DataType,HandleType,CoordinateType>::distributed_transfer()
{
    // Copy data from the source into a Tpetra vector.
    Tpetra::Vector<DataType> data_source_vector( 
	d_source_map, d_source->send_data(d_field_name));

    // Setup a Tpetra vector for the target. All entries set to zero.
    Tpetra::Vector<DataType> data_target_vector(d_target_map);

    // Transfer the data by exporting it from the source to the target.
    data_target_vector.doExport(data_source_vector, *d_export, Tpetra::INSERT);

    // Copy it to the target.
    data_target_vector.get1dCopy( d_target->receive_data(d_field_name) );
}

//---------------------------------------------------------------------------//

} // end namespace Coupler

#endif // core_Coupler_Data_Field_t_hpp

//---------------------------------------------------------------------------//
//                 end of Coupler_Data_Field.t.hpp
//---------------------------------------------------------------------------//
