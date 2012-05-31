//---------------------------------------------------------------------------//
/*!
 * \file DTK_Mesh.cpp
 * \author Stuart R. Slattery
 * \brief Concrete mesh definition for DTK algorithms.
 */
//---------------------------------------------------------------------------//

#include "DTK_Mesh.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Mesh::Mesh( const RCP_Moab& moab, 
	    const moab::Range& vertices,
	    const moab::Range& elements )
    : d_moab( moab )
    , d_vertices( vertices )
    , d_elements( elements )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Mesh::~Mesh()
{ /* ... */ }

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Mesh.cpp
//---------------------------------------------------------------------------//

