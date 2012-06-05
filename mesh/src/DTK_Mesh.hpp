//---------------------------------------------------------------------------//
/*!
 * \file DTK_Mesh.hpp
 * \author Stuart R. Slattery
 * \brief Concrete mesh declaration for DTK algorithms.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESH_HPP
#define DTK_MESH_HPP

#include <map>

#include <MBInterface.hpp>
#include <MBRange.hpp>

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

template<typename ElementHandle>
class Mesh
{

  public:

    //@{
    //! Typedefs.
    typedef ElementHandle                               element_handle_type;
    typedef Teuchos::RCP<moab::Interface>               RCP_Moab;
    typedef std::map<moab::EntityHandle,ElementHandle>  HandleMap;
    //@}

    // Constructor.
    Mesh( const RCP_Moab& moab, 
	  const moab::Range& elements,
	  const HandleMap& handle_map );

    // Destructor.
    ~Mesh();

    //! Get the Moab interface.
    const RCP_Moab& getMoab() const
    { return d_moab; }

    //! Get the mesh elements.
    const moab::Range& getElements() const
    { return d_elements; }

    // Given a moab element handle return the corresponding native element
    // handle.
    ElementHandle getNativeHandle( moab::EntityHandle moab_handle )
    { return d_handle_map[ moab_handle ]; }

  private:

    //! Moab interface implementation.
    RCP_Moab d_moab;

    //! Mesh elements.
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

// Create a mesh from an object the implements MeshTraits.
template<typename MeshObject>
Teuchos::RCP< Mesh<typename MeshObject::handle_type> > 
createMesh( const MeshObject& mesh_object );

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_Mesh_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_MESH_HPP

//---------------------------------------------------------------------------//
// end DTK_Mesh.hpp
//---------------------------------------------------------------------------//

