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
    // Generate the map for the data target.
    RCP_Tpetra_Map target_map 
	= Teuchos::rcp(
	    new Tpetra_Map_t(
		-1,
		transfer_data_field->target()->set_points( transfer_data_field->name() ),
		0,
		comm_global) );

    // Communicate points to the data source.
    int local_size 
	= transfer_data_field->target()->set_points( transfer_data_field->name() ).size();
    int local_max = 0;
    Teuchos::reduceAll<OrdinalType,int>(*comm_global,
					Teuchos::REDUCE_MAX, 
					int(1), 
					&local_size, 
					&local_max);
}

//---------------------------------------------------------------------------//

} // end namespace coupler

#endif // core_Mapper_t_hh

//---------------------------------------------------------------------------//
//                 end of Mapper.t.hh
//---------------------------------------------------------------------------//
