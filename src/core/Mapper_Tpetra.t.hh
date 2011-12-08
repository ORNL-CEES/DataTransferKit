//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Mapper.t.hh
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:31:05 2011
 * \brief  Mapper class template member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef core_Mapper_t_hh
#define core_Mapper_t_hh

#include <algorithm>

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR and DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class DataType, class HandleType, class CoordinateType>
Mapper<DataType,HandleType,CoordinateType>::Mapper()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class DataType, class HandleType, class CoordinateType>
Mapper<DataType,HandleType,CoordinateType>::~Mapper()
{ /* ... */ }


//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Map the field from the source onto the target.
 * \param comm_global The global communicator that the mapper will operate
 * on.
 * \param transfer_data_field The Transfer_Data_Field that will be
 * mapped.
 */
template<class DataType, class HandleType, class CoordinateType>
void Mapper<DataType,HandleType,CoordinateType>::map(
    const RCP_Communicator comm_global,
    RCP_Transfer_Data_Field transfer_data_field)
{
    // Get the local list of handles. These are the global indices for the
    // Tpetra map.
    const Teuchos::ArrayView<PointType> target_points = 
	transfer_data_field->target()->set_points( transfer_data_field->name() );
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
    RCP_Tpetra_Map data_target_map = 
	Tpetra::createNonContigMap<HandleType>( target_handles_view, comm_global);

    // Data target communicate points to the data source.
    int local_size 
	= transfer_data_field->target()->set_points( 
	    transfer_data_field->name() ).size();
    int local_max = 0;
    Teuchos::reduceAll<OrdinalType,int>(*comm_global,
					Teuchos::REDUCE_MAX, 
					int(1), 
					&local_size, 
					&local_max);
    
    // The data source finds the points in its domain.
    std::vector<HandleType> source_handles;


    
    // Generate the map for the data source.
    const Teuchos::ArrayView<const HandleType> source_handles_view(source_handles);
    RCP_Tpetra_Map data_source_map = 
	Tpetra::createNonContigMap<HandleType>( source_handles_view, comm_global);

    // Set the mapping in the field.
    transfer_data_field->set_mapping(data_source_map, data_target_map);
}

//---------------------------------------------------------------------------//

} // end namespace coupler

#endif // core_Mapper_t_hh

//---------------------------------------------------------------------------//
//                 end of Mapper.t.hh
//---------------------------------------------------------------------------//
