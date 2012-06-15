//---------------------------------------------------------------------------//
/*!
 * \file DTK_TopologyTools.hpp
 * \author Stuart R. Slattery
 * \brief TopologyTools declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TOPOLOGYTOOLS_HPP
#define DTK_TOPOLOGYTOOLS_HPP

#include <MBInterface.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class TopologyTools
 * \brief Tools based on concrete mesh topologies.
 */ 
//---------------------------------------------------------------------------//
class TopologyTools
{
  public:

    //! Constructor.
    TopologyTools()
    { /* ... */ }

    //! Destructor.
    ~TopologyTools()
    { /* ... */ }

    // Get the number of linear nodes for a particular iMesh topology.
    static int numLinearNodes( const moab::EntityType element_topology );

    // Point in element query.
    static bool pointInElement( Teuchos::Array<double>& coords,
				const moab::EntityHandle element,
				const Teuchos::RCP<moab::Interface>& moab );
};

} // end namepsace DataTransferKit

#endif // end DTK_TOPOLOGYTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_TopologyTools.hpp
//---------------------------------------------------------------------------//
