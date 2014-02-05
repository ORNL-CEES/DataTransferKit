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
 * \file DTK_ElementTree.hpp
 * \author Stuart R. Slattery
 * \brief ElementTree declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ELEMENTTREE_HPP
#define DTK_ELEMENTTREE_HPP

#include <boost/tr1/unordered_map.hpp>

#include "DTK_MeshManager.hpp"
#include "DTK_IntrepidCell.hpp"
#include "DTK_CloudSearch.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class ElementTree 
 * \brief A ElementTree data structure for local mesh searching in the
 * rendezvous decomposition.
 */
//---------------------------------------------------------------------------//
template<class Mesh>
class ElementTree
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                                mesh_type;
    typedef MeshTraits<Mesh>                                    MT;
    typedef typename MT::global_ordinal_type                    GlobalOrdinal;
    //@}

    // Constructor.
    ElementTree( const Teuchos::RCP<MeshManager<Mesh> >& mesh );

    // Destructor.
    ~ElementTree();

    // Build the kD-tree.
    void build();

    // Find a point in the tree.
    bool findPoint( const Teuchos::ArrayView<const double>& coords,
		    GlobalOrdinal& element,
		    double tolerance = 10*Teuchos::ScalarTraits<double>::eps() );

  private:

    // Mesh manager.
    Teuchos::RCP<MeshManager<Mesh> > d_mesh;

    // Local mesh element centroids.
    Teuchos::Array<double> d_element_centroids;

    // Vertex global-to-local indexers for each block.
    Teuchos::Array<boost::tr1::unordered_map<GlobalOrdinal,int> > d_vertex_g2l;

    // Discretization cell for each block.
    Teuchos::Array<IntrepidCell<Intrepid::FieldContainer<double> > > d_dcells;

    // Array of cumulative number elements in each block.
    Teuchos::Array<GlobalOrdinal> d_cumulative_elements;

    // Cloud Search.
    Teuchos::RCP<SearchTree> d_tree;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_ElementTree_impl.hpp"

//---------------------------------------------------------------------------//

#endif // DTK_ELEMENTTREE_HPP

//---------------------------------------------------------------------------//
// end ElementTree.hpp
//---------------------------------------------------------------------------//

