//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
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
#include <Teuchos_ArrayRCP.hpp>
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
    typedef Mesh                                mesh_type;
    typedef MeshTraits<Mesh>                    MT;
    typedef typename MT::global_ordinal_type    GlobalOrdinal;
    typedef Teuchos::Comm<int>                  CommType;
    typedef Teuchos::RCP<const CommType>        RCP_Comm;
    //@}

    //! Constructor.
    MeshTools()
    { /* ... */ }

    //! Destructor.
    ~MeshTools()
    { /* ... */ }


    //@{
    //! Bounds-checking mesh data access methods.
    // Get a view of the of the mesh nodes. The ArrayRCP object will not
    // manage the memory. 
    static Teuchos::ArrayRCP<const GlobalOrdinal> 
    nodesView( const Mesh& mesh );

    // Get a non-const view of the of the mesh nodes. The ArrayRCP object will
    // not manage the memory.
    static Teuchos::ArrayRCP<GlobalOrdinal> 
    nodesNonConstView( const Mesh& mesh );

    // Get a view of the of the mesh coordinates. The ArrayRCP object will not
    // manage the memory. 
    static Teuchos::ArrayRCP<const double> coordsView( const Mesh& mesh );

    // Get a non-const view of the of the mesh coordinates. The ArrayRCP
    // object will not manage the memory.
    static Teuchos::ArrayRCP<double> coordsNonConstView( const Mesh& mesh );

    // Get a view of the of the mesh elements. The ArrayRCP object will not
    // manage the memory. 
    static Teuchos::ArrayRCP<const GlobalOrdinal> 
    elementsView( const Mesh& mesh );

    // Get a non-const view of the of the mesh elements. The ArrayRCP object
    // will not manage the memory.
    static Teuchos::ArrayRCP<GlobalOrdinal> 
    elementsNonConstView( const Mesh& mesh );

    // Get a view of the of the mesh connectivity. The ArrayRCP object will not
    // manage the memory. 
    static Teuchos::ArrayRCP<const GlobalOrdinal> 
    connectivityView( const Mesh& mesh );

    // Get a non-const view of the of the mesh connectivity. The ArrayRCP
    // object will not manage the memory.
    static Teuchos::ArrayRCP<GlobalOrdinal> 
    connectivityNonConstView( const Mesh& mesh );

    // Get a view of the of the mesh connectivity permutation list. The
    // ArrayRCP object will not manage the memory. 
    static Teuchos::ArrayRCP<const std::size_t> 
    permutationView( const Mesh& mesh );

    // Get a non-const view of the of the mesh connectivity permutation
    // list. The ArrayRCP object will not manage the memory.
    static Teuchos::ArrayRCP<std::size_t> 
    permutationNonConstView( const Mesh& mesh );
    //@}


    //@{
    //! Bounding box methods.
    // Get the local bounding box for a mesh.
    static BoundingBox localBoundingBox( const Mesh& mesh );

    // Get the global bounding box for a mesh.
    static BoundingBox globalBoundingBox( const Mesh& mesh, 
					  const RCP_Comm& comm );
    //@}
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
