//---------------------------------------------------------------------------//
/*!
 * \file DTK_CellTopologyFactory.hpp
 * \author Stuart R. Slattery
 * \brief Factory method declaration for cell topology data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CELLTOPOLOGYFACTORY_HPP
#define DTK_CELLTOPOLOGYFACTORY_HPP

#include <MBInterface.hpp>

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{

class CellTopologyFactory
{
  public:

    // Consructor.
    CellTopologyFactory();

    // Destructor.
    ~CellTopologyFactory();

    // Factory method.
    Teuchos::RCP<shards::CellTopology> create( 
	const moab::EntityType element_topology, const int num_element_nodes );
};

} // end namespace DataTransferKit

#endif // end DTK_CELLTOPOLOGYFACTORY_HPP

//---------------------------------------------------------------------------//
// end DTK_CellTopologyFactory.hpp
//---------------------------------------------------------------------------//
