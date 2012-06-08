//---------------------------------------------------------------------------//
/*!
 * \file DTK_InverseComm_def.hpp
 * \author Stuart R. Slattery
 * \brief Inverse communication operator definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INVERSECOMM_DEF_HPP
#define DTK_INVERSECOMM_DEF_HPP

#include <vector>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Tpetra_Distributor.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename Ordinal>
InverseComm<Ordinal>::InverseComm( const RCP_Comm& comm )
: d_comm( comm )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*! 
 * \brief Destructor.
 */
template<typename Ordinal>
InverseComm<Ordinal>::~InverseComm()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Given a multimap of target procs and ordinals, compute the
 * corresponding Tpetra map. 
 */
template<typename Ordinal>
InverseComm<Ordinal>::RCP_TpetraMap 
InverseComm<Ordinal>::createImportMap( const std::multimap<int,Ordinal>& data )
{
    // Unroll the destinations for messages to be sent from this proccess.
    std::vector<int> data_dest;
    std::multimap<int,Ordinal>::const_iterator data_iterator;
    for ( data_iterator = data.begin(); data_iterator != data.end(); 
	  ++data_iterator )
    {
	for ( int n = 0; n < data.count( data_iterator->first ); ++n )
	{
	    data_dest.push_back( data_iterator->first );
	}
    }

    // Create a distributor.
    Tpetra::Distributor distributor( d_comm );
    distributor.createFromSends( Teuchos::ArrayView<const int>( data_dest ) );

    // Build the map from the distributor. First we'll post receives, then
    // we'll do the sends.
    Teuchos::ArrayView<const int> source_procs = distributor->getImagesFrom();
    Teuchos::ArrayView<const int> num_data = distributor->getLengthsFrom();
    Teuchos::ArrayView<const int>::const_iterator source_procs_iterator;
    Teuchos::ArrayView<const int>::const_iterator num_data_iterator;
    for ( source_procs_iterator = source_procs.begin(),
	      num_data_iterator = num_data.begin();
	  source_procs_iterator != source_procs.end();
	  ++source_procs_iterator, ++num_data_iterator )
    {
	
    }
    
    Teuchos::ArrayView<Ordinal> map_ids_view( map_ids );
    RCP_TpetraMap tpetra_map = 
	Tpetra::createNonContigMap<Ordinal>( map_ids, d_comm );

    testPostcondition( tpetra_map != Teuchos::null,
		       "Error creating communication map." );

    return tpetra_map;
}

} // end namespace DataTransferKit

#endif // end DTK_INVERSECOMM_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_InverseComm_def.hpp
//---------------------------------------------------------------------------//

