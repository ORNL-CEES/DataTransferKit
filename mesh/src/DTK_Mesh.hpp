//---------------------------------------------------------------------------//
/*!
 * \file DTK_Mesh.hpp
 * \author Stuart R. Slattery
 * \brief Concrete mesh declaration for DTK algorithms.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESH_HPP
#define DTK_MESH_HPP

#include <DTK_DataSource.hpp>

#include <MBInterface.hpp>
#include <MBRange.hpp>

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

class Mesh
{

  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<moab::Interface> RCP_Moab;
    //@}

    // Constructor.
    Mesh( const RCP_Moab& moab, 
	  const moab::Range& vertices,
	  const moab::Range& elements );

    // Destructor.
    ~Mesh();

    //! Get the Moab interface.
    RCP_Moab getMoab() const
    { return d_moab; }

    //! Get the mesh vertices.
    const moab::Range& getVertices() const
    { return d_vertices; }

    //! Get the mesh elements.
    const moab::Range& getElements() const
    { return d_elements; }

  private:

    //! Moab interface implementation.
    RCP_Moab d_moab;

    //! Mesh vertices.
    moab::Range d_vertices;

    //! Mesh elements.
    moab::Range d_elements;
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
// Non-member creation functions.
// Create a mesh from a DataSource.
template<typename NodeField, typename ElementField, typename DataField>
Teuchos::RCP<Mesh> createMeshFromDataSource( 
    const Teuchos::RCP< 
    DataSource<NodeField,ElementField,DataField> >& data_source );

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

