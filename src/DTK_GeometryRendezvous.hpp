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
 * \file DTK_GeometryRendezvous.hpp
 * \author Stuart R. Slattery
 * \brief Rendezvous decomposition for geometry declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYRENDEZVOUS_HPP
#define DTK_GEOMETRYRENDEZVOUS_HPP

#include <boost/tr1/unordered_map.hpp>

#include "DTK_GeometryTraits.hpp"
#include "DTK_GeometryManager.hpp"
#include "DTK_Partitioner.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Directory.hpp>


namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 \class GeometryRendezvous
 \brief Rendezvous decomposition for parallel geometry searching.

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
template<class Geometry, class GlobalOrdinal>
class GeometryRendezvous
{
  public:

    //@{
    //! Typedefs.
    typedef Geometry                                    geometry_type;
    typedef GlobalOrdinal                               global_ordinal_type;
    typedef GeometryTraits<Geometry>                    GT;
    typedef GeometryManager<Geometry,GlobalOrdinal>     GeometryManagerType;
    typedef Teuchos::RCP<GeometryManagerType>           RCP_GeometryManager;
    typedef Teuchos::RCP<Partitioner>                   RCP_Partitioner;
    typedef Teuchos::Comm<int>                          CommType;
    typedef Teuchos::RCP<const CommType>                RCP_Comm;
    //@}

    // Constructor.
    GeometryRendezvous( const RCP_Comm& comm, const int dimension,
			const BoundingBox& global_box );

    // Destructor.
    ~GeometryRendezvous();

    // Build the rendezvous decomposition.
    void build( const RCP_GeometryManager& geometry_manager );

    // Get the rendezvous destination processes for a blocked list of vertex
    // coordinates that are in the primary decomposition.
    Teuchos::Array<int> 
    procsContainingPoints( const Teuchos::ArrayRCP<double>& coords ) const;

    // Get the rendezvous destination processes for a set of bounding boxes.
    Teuchos::Array<Teuchos::Array<int> >
    procsContainingBoxes( const Teuchos::Array<BoundingBox>& boxes ) const;

    // Get the native geometry ids in the rendezvous decomposition and their
    // source decomposition procs containing a blocked list of coordinates
    // also in the rendezvous decomposition.
    void geometryContainingPoints( 
	const Teuchos::ArrayRCP<double>& coords,
	Teuchos::Array<GlobalOrdinal>& gids,
	Teuchos::Array<int>& geometry_src_procs,
	const double geometric_tolerance ) const;

    //! Get the bounding box over which the rendezvous decomposition was
    //! generated.
    const BoundingBox& getBox() const
    { return d_global_box; }

    //! For a list of geometry gids in the rendezvous decomposition, get their
    //! source procs.
    Teuchos::Array<int> geometrySourceProcs( 
	const Teuchos::Array<GlobalOrdinal>& gids );

  private:

    // Extract the geometry that is in a bounding box.
    void getGeometryInBox( const RCP_GeometryManager& geometry_manager );

    // Send the geometry to the rendezvous decomposition.    
    void sendGeometryToRendezvous( const RCP_GeometryManager& geometry_manager );

  private:

    // Global communicator over which to perform the rendezvous.
    RCP_Comm d_comm;

    // The dimension of the rendezvous.
    int d_dimension;

    // Bounding box in which to perform the rendezvous.
    BoundingBox d_global_box;

    // Rendezvous partitioning.
    RCP_Partitioner d_partitioner;

    // Rendezvous geometry gid to source proc map.
    std::tr1::unordered_map<GlobalOrdinal,int> d_geometry_src_procs_map;

    // Rendezvous on-process geometry.
    Teuchos::Array<Geometry> d_rendezvous_geometry;

    // Rendezvous on-process geometry gids.
    Teuchos::Array<GlobalOrdinal> d_rendezvous_gids;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_GeometryRendezvous_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_GEOMETRYRENDEZVOUS_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryRendezvous.hpp
//---------------------------------------------------------------------------//

