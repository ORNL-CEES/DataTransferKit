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
 * \file DTK_KDTree.hpp
 * \author Stuart R. Slattery
 * \brief KDTree declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_KDTREE_HPP
#define DTK_KDTREE_HPP

#include "DTK_RendezvousMesh.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <MBAdaptiveKDTree.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class KDTree
 * \brief A KDTree data structure for local mesh searching in the rendezvous
 * decomposition.

 A kD-tree provides n*log(n) construction time and space complexity with
 log(n) search time complexity.

 */
//---------------------------------------------------------------------------//
template<typename GlobalOrdinal>
class KDTree
{
  public:

    //@{
    //! Typedefs.
    typedef GlobalOrdinal                        global_ordinal_type;
    typedef RendezvousMesh<GlobalOrdinal>        RendezvousMeshType;
    typedef Teuchos::RCP<RendezvousMeshType>     RCP_RendezvousMesh;
    //@}

    // Constructor.
    KDTree( const RCP_RendezvousMesh& mesh, const int dim );

    // Destructor.
    ~KDTree();

    // Build the kD-tree.
    void build();

    // Find a point in the tree.
    bool findPoint( const Teuchos::Array<double>& coords,
		    GlobalOrdinal& element );

  private:

    // Find a point in a leaf.
    bool findPointInLeaf( const Teuchos::Array<double>& coords,
			  const moab::EntityHandle leaf,
			  moab::EntityHandle& element );

  private:

    // Moab Mesh.
    RCP_RendezvousMesh d_mesh;

    // Tree dimension.
    int d_dim;

    // Adaptive kD-tree.
    moab::AdaptiveKDTree d_tree;

    // Tree root.
    moab::EntityHandle d_root;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_KDTree_def.hpp"

#endif // DTK_KDTREE_HPP

//---------------------------------------------------------------------------//
// end KDTree.hpp
//---------------------------------------------------------------------------//

