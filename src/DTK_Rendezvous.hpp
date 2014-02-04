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
 * \file DTK_Rendezvous.hpp
 * \author Stuart R. Slattery
 * \brief Rendezous declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUS_HPP
#define DTK_RENDEZVOUS_HPP

#include <boost/tr1/unordered_map.hpp>

#include "DTK_MeshTraits.hpp"
#include "DTK_MeshManager.hpp"
#include "DTK_RendezvousMesh.hpp"
#include "DTK_MeshContainer.hpp"
#include "DTK_KDTree.hpp"
#include "DTK_Partitioner.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Directory.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 \class Rendezvous
 \brief Rendezvous decomposition for parallel mesh searching.

 Relating two non-conformal meshes will ultimately require some type of
 evaluation algorithm to apply the data from one geometry to another. To drive
 these evaluation algorithms, the target objects to which this data will be
 applied must be located within the the source geometry. In a serial
 formulation, efficient search structures that offer logarithmic asymptotic
 time complexity are available to perform this operation. However, in a
 parallel formulation, if these two geometries are arbitrarily decomposed,
 geometric alignment is not likely and a certain degree of communication will
 be required. A geometric rendezvous manipulates the source and target
 geometries such that all geometric operations have a local formulation.

 A geometry that is associated with the providing data through function
 evaluations will be referred to as the source geometry while the geometry
 that will be receiving the data will be referred to as the target
 geometry. The rendezvous decomposition has several properties. It is defined
 over a communicator that encapsulates the union of the communication spaces
 owned by the source and target geometries. It is defined inside of a global,
 axis-aligned bounding box that bounds the intersection of the source and
 target geometries. The decomposition is of the same dimension as the source
 and target geometries. A rendezvous decomposition cannot be generated with
 source and target geometries of different dimensions (e.g. a 3 dimensional
 source geometry and a 2 dimensional target geometry cannot be used to
 generate a rendezvous decomposition).
 */
//---------------------------------------------------------------------------//
template<class Mesh>
class Rendezvous
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                        mesh_type;
    typedef MeshTraits<Mesh>                            MT;
    typedef typename MT::global_ordinal_type            GlobalOrdinal;
    typedef Teuchos::RCP< MeshManager<Mesh> >           RCP_MeshManager;
    typedef typename MeshManager<Mesh>::BlockIterator   BlockIterator;
    typedef MeshContainer<GlobalOrdinal>                MeshContainerType;
    typedef RendezvousMesh<GlobalOrdinal>               RendezvousMeshType;
    typedef Teuchos::RCP<RendezvousMeshType>            RCP_RendezvousMesh;
    typedef KDTree<GlobalOrdinal>                       KDTreeType;
    typedef Teuchos::RCP<KDTreeType>                    RCP_KDTree;
    typedef Teuchos::RCP<Partitioner>                   RCP_Partitioner;
    typedef Teuchos::Comm<int>                          CommType;
    typedef Teuchos::RCP<const CommType>                RCP_Comm;
    typedef Tpetra::Map<int,GlobalOrdinal>              TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>               RCP_TpetraMap;
    //@}

    // Constructor.
    Rendezvous( const RCP_Comm& comm, const int dimension,
		const BoundingBox& global_box );

    // Destructor.
    ~Rendezvous();

    // Build the rendezvous decomposition.
    void build( const RCP_MeshManager& mesh_manager );

    // Get the rendezvous destination processes for a blocked list of vertex
    // coordinates that are in the primary decomposition.
    Teuchos::Array<int> 
    procsContainingPoints( const Teuchos::ArrayRCP<double>& coords ) const;

    // Get the rendezvous destination processes for a set of bounding boxes.
    Teuchos::Array<Teuchos::Array<int> >
    procsContainingBoxes( const Teuchos::Array<BoundingBox>& boxes ) const;

    // Get the native mesh elements in the rendezvous decomposition and their
    // source decomposition procs containing a blocked list of coordinates
    // also in the rendezvous decomposition.
    void elementsContainingPoints( 
	const Teuchos::ArrayRCP<double>& coords,
	Teuchos::Array<GlobalOrdinal>& elements,
	Teuchos::Array<int>& element_src_procs,
	double tolerance = 10*Teuchos::ScalarTraits<double>::eps() ) const;

    // Get the native elements in the rendezvous decomposition that are in
    // each bounding box in a list.
    void elementsInBoxes(
	const Teuchos::Array<BoundingBox>& boxes,
	Teuchos::Array<Teuchos::Array<GlobalOrdinal> >& elements ) const;

    // Get the native elements in the rendezvous decomposition that are in
    // each geometry in a list.
    template<class Geometry>
    void elementsInGeometry(
	const Teuchos::Array<Geometry>& geometry,
	Teuchos::Array<Teuchos::Array<GlobalOrdinal> >& elements,
	const double tolerance,	bool all_vertices_for_inclusion ) const;

    //! Get the rendezvous mesh.
    const RCP_RendezvousMesh& getMesh() const
    { return d_rendezvous_mesh; }

    //! Get the bounding box over which the rendezvous decomposition was
    //! generated.
    const BoundingBox& getBox() const
    { return d_global_box; }

    //! For a list of elements in the rendezvous decomposition, get their
    //! source procs.
    Teuchos::Array<int> elementSourceProcs( 
	const Teuchos::Array<GlobalOrdinal>& elements );

  private:

    // Extract the mesh block vertices and elements that are in a bounding box.
    void getMeshInBox( const RCP_MeshManager& mesh_manager );

    // Send the mesh to the rendezvous decomposition and build the concrete
    // mesh blocks.
    MeshManager<MeshContainerType> 
    sendMeshToRendezvous( const RCP_MeshManager& mesh_manager );

    // Setup the import communication patterns.
    void setupImportCommunication( 
	const Teuchos::RCP<Mesh>& mesh,
	const Teuchos::ArrayView<short int>& elements_in_box,
	Teuchos::Array<GlobalOrdinal>& rendezvous_vertices,
	Teuchos::Array<GlobalOrdinal>& rendezvous_elements );

  private:

    // Global communicator over which to perform the rendezvous.
    RCP_Comm d_comm;

    // The dimension of the rendezvous.
    int d_dimension;

    // Bounding box in which to perform the rendezvous.
    BoundingBox d_global_box;

    // Rendezvous partitioning.
    RCP_Partitioner d_partitioner;

    // Rendezvous mesh element to source proc map.
    std::tr1::unordered_map<GlobalOrdinal,int> d_element_src_procs_map;

    // Rendezvous on-process mesh.
    RCP_RendezvousMesh d_rendezvous_mesh;

    // Rendezvous on-process kD-tree.
    RCP_KDTree d_kdtree;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_Rendezvous_def.hpp"

#endif // end DTK_RENDEZVOUS_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous.hpp
//---------------------------------------------------------------------------//

