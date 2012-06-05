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
#include <Teuchos_ENull.hpp>

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
    // Extract the local source nodes and elements that are within the
    // provided bounding box. 
    std::vector<NodeType> local_nodes;
    std::vector<ElementType> local_elements;
    getMeshInBox( source_nodes, source_elements, local_nodes, local_elements );
        
    // Construct the rendezvous decomposition of the mesh with RCB.
    d_rcb = 
	Teuchos::rcp( new RCB<SourceMeshNodeField>( local_nodes, d_comm ) );

    testPostcondition( d_rcb != Teuchos::null,
		       "Error creating RCB decomposition." );

    d_rcb->partition();

    // Send the mesh description to the rendezvous decomposition.

    // Clear the local copy of the source mesh.
    local_nodes.clear();
    local_elements.clear();

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
/*!
 * \brief Get the source mesh elements containing a list of coordinates.
 */
template<typename SourceMeshNodeField, typename SourceMeshElementField>
Rendezvous<SourceMeshNodeField,SourceMeshElementField>::SourceElementHandleVec
Rendezvous<SourceMeshNodeField,SourceMeshElementField>::getElementHandles( 
    const std::vector<double>& coords ) const
{
    testPrecondition( coords.size() % 3 == 0, 
		      "Three dimensional coordinates not provided." );
    
    int num_points = coords.size() / 3;
    SourceElementHandleVec element_handles;
    for ( int i = 0; i < num_points; ++i )
    {
	double point[3] = { coords[3*i], coords[3*i+1], coords[3*i+2] };
	source_element_handle_type handle = d_kdtree->findPoint( point );
	element_handles.push_back( handle );
    }

    testPostcondition( (int) element_handles.size() == num_points,
		       "Error getting source mesh elements." );

    return element_handles;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extract the mesh nodes and elements that are in the bounding box.
 */
template<typename SourceMeshNodeField, typename SourceMeshElementField>
void
Rendezvous<SourceMeshNodeField,SourceMeshElementField>::getMeshInBox(
    const SourceMeshNodeField& source_nodes,
    const SourceMeshElementField& source_elements,
    std::vector<NodeType>& local_nodes,
    std::vector<ElementType>& local_elements )
{

}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//

