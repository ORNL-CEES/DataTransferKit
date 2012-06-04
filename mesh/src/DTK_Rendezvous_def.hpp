//---------------------------------------------------------------------------//
/*!
 * \file DTK_Rendezvous_def.hpp
 * \author Stuart R. Slattery
 * \brief Rendezvous definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUS_DEF_HPP
#define DTK_RENDEZVOUS_DEF_HPP

#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename SourceMeshNodeField, typename SourceMeshElementField>
Rendezvous<SourceMeshNodeField,SourceMeshElementField>::Rendezvous(
    const RCP_Comm& global_comm, const BoundingBox& global_box )
    : d_global_comm( global_comm )
    , d_global_box( global_box )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename SourceMeshNodeField, typename SourceMeshElementField>
Rendezvous<SourceMeshNodeField,SourceMeshElementField>::~Rendezvous()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the rendezvous decomposition.
 */
template<typename SourceMeshNodeField, typename SourceMeshElementField>
Rendezvous<SourceMeshNodeField,SourceMeshElementField>::build(
    const SourceMeshNodeField& source_nodes, 
    const SourceMeshElementField& source_elements )
{
    // Construct the rendezvous decomposition of the mesh with RCB.
    d_rcb = 
	Teuchos::rcp( new RCB<SourceMeshNodeField>( source_nodes, d_comm ) );
    d_rcb->partition();

    // Send the mesh description to the rendezvous decomposition.

    // Build the concrete mesh database in the rendezvous decomposition.
    
    // Create a kD-tree in the rendezvous decomposition.
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get the rendezvous process for a list of points.
 */
template<typename SourceMeshNodeField, typename SourceMeshElementField>
std::vector<int> 
Rendezvous<SourceMeshNodeField,SourceMeshElementField>::getRendezvousProc(
    const std::vector<double>& coords ) const
{
    testPrecondition( coords.size() % 3 == 0, 
		      "Three dimensional coordinates not provided." );
    
    int num_points = coords.size() / 3;
    std::vector<int> destination_procs;
    for ( int i = 0; i < num_points; ++i )
    {
	double point[3] = { coords[3*i], coords[3*i+1], coords[3*i+2] };
	int rendezvous_proc = d_rcb->getDestinationProc( point );
	destination_procs.push_back( rendezvous_proc );
    }

    testPostcondition( (int) destination_procs.size() == num_points,
		       "Error getting destination processes." );

    return destination_procs;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//

