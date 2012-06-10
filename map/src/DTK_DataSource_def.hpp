//---------------------------------------------------------------------------//
/*!
 * \file DTK_DataSource_def.hpp
 * \author Stuart R. Slattery
 * \brief DataSource defintion.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATASOURCE_DEF_HPP
#define DTK_DATASOURCE_DEF_HPP

#include <DTK_Exception.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor.
 * \param mpi_comm Raw MPI communicator over which this data source is
 * defined. 
 * \param mesh The mesh for this data source. This object must implement
 * MeshTraits.
 */
template<class Mesh, class DataField>
DataSource<Mesh,DataField>::DataSource( const Mesh& mesh,
					const MPI_Comm& mpi_comm )
    : d_mesh( mesh )
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
template<class Mesh, class DataField>
DataSource<Mesh,DataField>::~DataSource()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Register a source field evaluator.
 * \param field_name The string key for this field.
 * \param field_evaluator The field evaluator for this field.
 */
template<class Mesh, class DataField>
void DataSource<Mesh,DataField>::registerFieldEvaluator( 
    const std::string& name field_name, const FieldEvaluator& field_evaluator )
{
    // Create an integer id for this field.
    std::size_t field_id = d_name_map.size();

    // Add it to the maps.
    d_name_map[ field_name ] = field_id;
    d_eval_map[ field_id ] = field_evaluator;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_DATASOURCE_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_DataSource_def.hpp
//---------------------------------------------------------------------------//
