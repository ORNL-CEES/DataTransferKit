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

namespace coupler
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
	map();
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
	distributed_transfer();
    }
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Generate topology map for this field.
 */
template<class DataType, class HandleType, class CoordinateType>
void Data_Field<DataType,HandleType,CoordinateType>::map()
{
    // Get the local list of handles. These are the global indices for the
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
    const Teuchos::ArrayView<const HandleType> target_handles_view(target_handles);
    d_target_map = 
	Tpetra::createNonContigMap<HandleType>( target_handles_view, d_comm);

    // Data target communicate points to the data source.
    int local_size 
	= d_target->set_points( d_field_name ).size();
    int local_max = 0;
    Teuchos::reduceAll<OrdinalType,int>(*d_comm,
					Teuchos::REDUCE_MAX, 
					int(1), 
					&local_size, 
					&local_max);
    
    // The data source finds the points in its domain.
    std::vector<HandleType> source_handles;


    
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

}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform distributed transfer.
 */
template<class DataType, class HandleType, class CoordinateType>
void Data_Field<DataType,HandleType,CoordinateType>::distributed_transfer()
{

}

//---------------------------------------------------------------------------//

} // end namespace coupler

#endif // core_Coupler_Data_Field_t_hpp

//---------------------------------------------------------------------------//
//                 end of Coupler_Data_Field.t.hpp
//---------------------------------------------------------------------------//
