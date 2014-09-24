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
 * \file DTK_MeshContainer.hpp
 * \author Stuart R. Slattery
 * \brief A simple default mesh container with a mesh block
 * implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHCONTAINER_HPP
#define DTK_MESHCONTAINER_HPP

#include "DTK_MeshBlock.hpp"
#include "DTK_MeshTypes.hpp"

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class MeshContainer
 * \brief A default mesh block implementation.

 The MeshContainer is a default implementation of MeshBlock for clients. Its
 data mimics the structure of MeshBlock with constructor arguments expected
 to have the data layout of a MeshBlock object. Each block of mesh must be in
 a separate mesh container.
 */
//---------------------------------------------------------------------------//
class MeshContainer : public MeshBlock
{
  public:

    //@{
    //! Typedefs.
    typedef MeshBlock Base;
    //@}
    
    //! Default Constructor.
    MeshContainer()
    { /* ... */ }

    /*!
     * \brief Constructor.
     *
     * \param vertex_dim The dimension of the vertices in the mesh.
     *
     * \param vertices The vertex global ordinals in the mesh block.
     *
     * \param coordinates The vertex coordinates in the mesh block. These must
     * be blocked by dimension.
     *
     * \param element_topology The DTK_ElementTopology of the elements in the
     * mesh block.
     *
     * \param vertices_per_element The number of vertices used to construct
     * each element in the mesh block.
     *
     * \param elements The element global ordinals in the mesh block.
     *
     * \param connectivity The vertex global ordinals that construct the
     * elements in the mesh block. The connectivity values must be blocked by
     * client canonical vertex ordering.
     *
     * \param permutation_list The permuation list describing the difference
     * between DTK canonical vertex ordering and client canonical vertex
     * ordering for the element topology in this mesh block. This list must be
     * the same size as vertices_per_element.
     */
    MeshContainer( 
	int dimension,
	const Teuchos::ArrayRCP<MeshId>& vertices,
	const Teuchos::ArrayRCP<double>& coords,
	DTK_ElementTopology element_topology,
	int vertices_per_element,
	const Teuchos::ArrayRCP<MeshId>& elements,
	const Teuchos::ArrayRCP<MeshId>& connectivity,
	const Teuchos::ArrayRCP<int>& permutation_list )
	: d_dimension( dimension )
	, d_vertices( vertices )
	, d_coords( coords )
	, d_element_topology( element_topology )
	, d_vertices_per_element( vertices_per_element )
	, d_elements( elements )
	, d_connectivity( connectivity )
	, d_permutation_list( permutation_list )
    { /* ... */ }

    //! Destructor.
    ~MeshContainer()
    { /* ... */ }

    /*!
     * \brief Return the spatial dimension of the block.
     */
    int dimension()
    { return d_dimension; }

    /*!
     * \brief Return a reference-counted pointer to the vertex ids.
     */
    Teuchos::ArrayRCP<MeshId> vertexIds()
    { return d_vertices; }

    /*!
     * \brief Return a reference-counted pointer to the vertex coordinates.
     */
    Teuchos::ArrayRCP<double> vertexCoordinates()
    { return d_coords; }

    /*! 
     * \brief Return the element topology for this mesh block
     * (DTK_ElementTopology enum).
     */
    DTK_ElementTopology elementTopology()
    { return d_element_topology; }

    /*! 
     * \brief Return the number of vertices that constructs an individual
     * element in this mesh block. All elements in the mesh must be
     * constructed with the same number of vertices.
     */
    int verticesPerElement()
    { return d_vertices_per_element; }

    /*! 
     * \brief Return a reference-counted pointer to the element ids.
     */
    Teuchos::ArrayRCP<MeshId> elementIds()
    { return d_elements; }

    /*! 
     * \brief Return a reference-counted pointer to the element connectivity.
     */
    Teuchos::ArrayRCP<MeshId> connectivity()
    { return d_connectivity; }

    /*! 
     * \brief Return a reference-counted pointer to the connectivity
     * permutation array.
     */
    Teuchos::ArrayRCP<int> permutation()
    { return d_permutation_list; }
    
  private:

    // Spatial dimension.
    int d_dimension;

    // Vertices.
    Teuchos::ArrayRCP<MeshId> d_vertices;

    // Coordinates.
    Teuchos::ArrayRCP<double> d_coords;

    // Element topology.
    DTK_ElementTopology d_element_topology;

    // Vertices per element.
    int d_vertices_per_element;

    // Elements.
    Teuchos::ArrayRCP<MeshId> d_elements;

    // Connectivity.
    Teuchos::ArrayRCP<MeshId> d_connectivity;

    // Permutation list.
    Teuchos::ArrayRCP<int> d_permutation_list;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESHCONTAINER_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshContainer.hpp
//---------------------------------------------------------------------------//
