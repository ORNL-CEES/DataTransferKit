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
#include <Teuchos_ArrayRCP.hpp>

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

    // Create a distributor. This data structure will provide this process
    // with knowledge of how many messages and of what sizes it will be
    // receiving and the source rank it will be receiving them from.
    Tpetra::Distributor distributor( d_comm );
    int num_imports = 
	distributor.createFromSends( Teuchos::arrayViewFromVector( data_dest ) );

    // Build the data packets.
    Teuchos::ArrayView<const int> source_procs = distributor->getImagesFrom();
    Teuchos::ArrayView<const int> num_data = distributor->getLengthsFrom();
    Teuchos::ArrayView<const int>::const_iterator source_procs_iterator;
    Teuchos::ArrayView<const int>::const_iterator num_data_iterator;
    std::vector<RCP_CommRequest> requests( num_data.size() );
    std::vector< Teuchos::ArrayRCP<Ordinal> > receive_packets;
    for ( source_procs_iterator = source_procs.begin(),
	      num_data_iterator = num_data.begin()
			  int i = 0;
	  source_procs_iterator != source_procs.end();
	  ++source_procs_iterator, ++num_data_iterator, ++i )
    {
	// Allocate a packet.
	receive_packets.push_back( 
	    Teuchos::ArrayRCP<Ordinal>( *num_data_iterator ) );
    }

    // Move them with the distributor.
    distributor.doPostsAndWaits();

    // Now we can make the map by unrolling the packets into an id vector.
    std::vector<Ordinal> map_ids;
    std::vector< Teuchos::ArrayRCP<Ordinal> >::const_iterator 
	receive_packet_iterator;
    Teuchos::ArrayRCP<Ordinal>::const_iterator ordinal_iterator;
    for ( receive_packet_iterator = receive_packets.begin();
	  receive_packet_iterator = receive_packets.end();
	  ++receive_packet_iterator )
    {
	for ( ordinal_iterator = receive_packet_iterator->begin();
	      ordinal_iterator = receive_packet_iterator->end();
	      ++ordinal_iterator )
	{
	    map_ids.push_back( *ordinal_iterator );
	}
    }
    receive_packets.clear();

    RCP_TpetraMap tpetra_map = Tpetra::createNonContigMap<Ordinal>( 
	Teuchos::arrayViewFromVector( map_ids ), d_comm );

    testPostcondition( tpetra_map != Teuchos::null,
		       "Error creating communication map." );

    return tpetra_map;
}

} // end namespace DataTransferKit

#endif // end DTK_INVERSECOMM_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_InverseComm_def.hpp
//---------------------------------------------------------------------------//

