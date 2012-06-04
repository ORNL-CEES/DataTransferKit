//---------------------------------------------------------------------------//
/*!
 * \file DTK_Mesh_def.hpp
 * \author Stuart R. Slattery
 * \brief Concrete mesh template definitions.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESH_DEF_HPP
#define DTK_MESH_DEF_HPP

#include <vector>
#include <cassert>

#include <DTK_Exception.hpp>
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include <MBCore.hpp>

#include <Teuchos_ENull.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename ElementHandle>
Mesh<ElementHandle>::Mesh( const RCP_Moab& moab, 
			   const moab::Range& elements,
			   const HandleMap& handle_map )
    : d_moab( moab )
    , d_elements( elements )
    , d_handle_map( handle_map )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename ElementHandle>
Mesh<ElementHandle>::~Mesh()
{ /* ... */ }

} // end namespace DataTransferKit

#endif // end DTK_MESH_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Mesh_def.hpp
//---------------------------------------------------------------------------//

