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

} // end namespace DataTransferKit

#endif // end DTK_DATATARGET_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_DataTarget_def.hpp
//---------------------------------------------------------------------------//
