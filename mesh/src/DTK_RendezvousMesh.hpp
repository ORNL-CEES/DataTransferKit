//---------------------------------------------------------------------------//
/*!
 * \file DTK_RendezvousMesh.hpp
 * \author Stuart R. Slattery
 * \brief Concrete rendezvous mesh declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUSMESH_HPP
#define DTK_RENDEZVOUSMESH_HPP

#include <map>

#include <DTK_MeshTraits.hpp>

#include <MBInterface.hpp>
#include <MBRange.hpp>

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

template<typename Handle>
class RendezvousMesh
{

  public:

    //@{
    //! Typedefs.
    typedef Handle                                    handle_type;
    typedef Teuchos::RCP<moab::Interface>             RCP_Moab;
    typedef std::map<moab::EntityHandle,handle_type>  HandleMap;
    //@}

    // Constructor.
    RendezvousMesh( const RCP_Moab& moab, 
		    const moab::Range& elements,
		    const HandleMap& handle_map );

    // Destructor.
    ~RendezvousMesh();

    //! Get the Moab interface.
    const RCP_Moab& getMoab() const
    { return d_moab; }

    //! Get the mesh elements.
    const moab::Range& getElements() const
    { return d_elements; }

    // Given a moab element handle return the corresponding native element
    // handle.
    handle_type getNativeHandle( moab::EntityHandle moab_handle )
    { return d_handle_map[ moab_handle ]; }

  private:

    //! Moab interface implementation.
    RCP_Moab d_moab;

    //! RendezvousMesh elements.
    moab::Range d_elements;

    //! Moab element handle to native element handle map.
    HandleMap d_handle_map;
};

//---------------------------------------------------------------------------//
//! DTK to MOAB topology translation table.
const moab::EntityType moab_topology_table[] =
{
    moab::MBVERTEX, // DTK_VERTEX
    moab::MBEDGE,   // DTK_LINE_SEGMENT
    moab::MBTRI,    // DTK_TRIANGLE
    moab::MBQUAD,   // DTK_QUADRILATERAL
    moab::MBTET,    // DTK_TETRAHEDRON
    moab::MBHEX     // DTK_HEXAHEDRON
};

//---------------------------------------------------------------------------//
// Non-member creation methods.
//---------------------------------------------------------------------------//

// Create a RendezvousMesh from an object that implements MeshTraits.
template<typename Mesh>
Teuchos::RCP< RendezvousMesh<typename MeshTraits<Mesh>::handle_type> > 
createRendezvousMesh( const Mesh& mesh );

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_RendezvousMesh_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_RENDEZVOUSESH_HPP

//---------------------------------------------------------------------------//
// end DTK_RendezvousMesh.hpp
//---------------------------------------------------------------------------//

