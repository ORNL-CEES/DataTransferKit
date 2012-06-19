//---------------------------------------------------------------------------//
/*!
 * \file DTK_MeshTools.hpp
 * \author Stuart R. Slattery
 * \brief MeshTools declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTOOLS_HPP
#define DTK_MESHTOOLS_HPP

#include "DTK_MeshTraits.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class MeshTools
 * \brief Tools for objects that have mesh traits.
 */ 
//---------------------------------------------------------------------------//
template<class Mesh>
class MeshTools
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                            mesh_type;
    typedef MeshTraits<Mesh>                MT;
    typedef Teuchos::Comm<int>              CommType;
    typedef Teuchos::RCP<const CommType>    RCP_Comm;
    //@}

    //! Constructor.
    MeshTools()
    { /* ... */ }

    //! Destructor.
    ~MeshTools()
    { /* ... */ }

    // Get the local bounding box for a mesh.
    static BoundingBox localBoundingBox( const Mesh& mesh );

    // Get the global bounding box for a mesh.
    static BoundingBox globalBoundingBox( const Mesh& mesh, 
					  const RCP_Comm& comm );
};

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_MeshTools_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_MESHTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshTools.hpp
//---------------------------------------------------------------------------//
