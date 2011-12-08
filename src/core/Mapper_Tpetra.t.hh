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
    const Teuchos::ArrayView<PointType> local_points = 
	transfer_data_field->target()->set_points( transfer_data_field->name() );
    typename Teuchos::ArrayView<PointType>::const_iterator point_it;

    std::vector<HandleType> local_handles(local_points.size());
    typename std::vector<HandleType>::iterator handle_it;

    for (handle_it = local_handles.begin(), point_it = local_points.begin(); 
	 handle_it != local_handles.end();
	 ++handle_it, ++point_it)
    {
	*handle_it = point_it->handle();
    }
    const Teuchos::ArrayView<const HandleType> local_handles_view(local_handles);

    // Generate the map for the data target.
    RCP_Tpetra_Map data_target_map = 
	Tpetra::createNonContigMap<HandleType>( local_handles_view, comm_global);

    // Communicate points to the data source.
    int local_size 
	= transfer_data_field->target()->set_points( transfer_data_field->name() ).size();
    int local_max = 0;
    Teuchos::reduceAll<OrdinalType,int>(*comm_global,
					Teuchos::REDUCE_MAX, 
					int(1), 
					&local_size, 
					&local_max);
    
    Teuchos::broadcast<OrdinalType,char>(*comm_global, 
					 comm_global->getRank(),
					 local_points);
					
    // The data source finds the points in its domain and builds its map.
					
}

//---------------------------------------------------------------------------//

} // end namespace coupler

#endif // core_Mapper_t_hh

//---------------------------------------------------------------------------//
//                 end of Mapper.t.hh
//---------------------------------------------------------------------------//
