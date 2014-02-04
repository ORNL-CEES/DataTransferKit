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
 * \file DTK_SharedDomainMap.hpp
 * \author Stuart R. Slattery
 * \brief Shared domain map declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SHAREDDOMAINMAP_HPP
#define DTK_SHAREDDOMAINMAP_HPP

#include <boost/tr1/unordered_map.hpp>

#include "DTK_FieldTraits.hpp"
#include "DTK_FieldEvaluator.hpp"
#include "DTK_FieldManager.hpp"
#include "DTK_MeshTraits.hpp"
#include "DTK_MeshManager.hpp"
#include "DTK_BoundingBox.hpp"
#include "DTK_CommIndexer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Directory.hpp>
#include <Tpetra_Export.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 \class SharedDomainMap
 \brief A map for shared domain problems.

 A parallel topology map is an operator, \f$M\f$, that defines the translation
 of a field, \f$F(s)\f$, from a source spatial domain, \f$\Omega_S\f$, to a
 field, \f$G(t)\f$, in the target spatial domain \f$\Omega_T\f$, such that
 \f$M: \mathcal{R}^D \rightarrow \mathcal{R}^D, \forall r \in [\Omega_S \cap
 \Omega_T]\f$, using both geometric and parallel operations.

 A shared domain problem is one in which the geometric domains of the source
 and target intersect over all dimensions of the problem. The shared domain
 map is a parallel topology map that has several properties. It is defined
 over a communicator that encapsulates the union of the communication spaces
 owned by the source and target geometries.  The map is of the same dimension
 as the source and target geometries. A shared domain map cannot be generated
 with source and target geometries of different dimensions (e.g. a 3
 dimensional source geometry and a 2 dimensional target geometry cannot be
 used to generate a shared domain map).

 There are several steps in the mapping where target objects may or may not be
 found in the source mesh. Both the RCB decomposition search and kD-tree search
 have the potential to return target objects that were not found in the source
 mesh. The source function will not be evaluated at these points as they are
 not in the domain \f$\Omega_S\f$, and therefore the evaluation operation will
 not be valid. However, a list of these points in the target decomposition may
 be generated for further processing by the client.

*/
//---------------------------------------------------------------------------//
template<class Mesh, class CoordinateField>
class SharedDomainMap
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                      mesh_type;
    typedef MeshTraits<Mesh>                          MT;
    typedef typename MT::global_ordinal_type          GlobalOrdinal;
    typedef MeshManager<Mesh>                         MeshManagerType;
    typedef Teuchos::RCP<MeshManagerType>             RCP_MeshManager;
    typedef typename MeshManagerType::BlockIterator   MeshBlockIterator;
    typedef CoordinateField                           coord_field_type;
    typedef FieldTraits<CoordinateField>              CFT;
    typedef typename CFT::size_type                   CoordOrdinal;
    typedef FieldManager<CoordinateField>             CoordFieldManagerType;
    typedef Teuchos::RCP<CoordFieldManagerType>       RCP_CoordFieldManager;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;
    typedef Tpetra::Map<int,GlobalOrdinal>            TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    typedef Tpetra::Export<int,GlobalOrdinal>         ExportType;
    typedef Teuchos::RCP<ExportType>                  RCP_TpetraExport;
    //!@}

    // Constructor.
    SharedDomainMap( const RCP_Comm& comm, const int dimension, 
		     bool store_missed_points = false );

    // Destructor.
    ~SharedDomainMap();

    // Generate the shared domain map.
    void setup( const RCP_MeshManager& source_mesh_manager, 
		const RCP_CoordFieldManager& target_coord_manager,
		double tolerance = 10*Teuchos::ScalarTraits<double>::eps() );

    // Apply the shared domain map by evaluating a function at target points
    // that were mapped.
    template<class SourceField, class TargetField>
    void apply( 
	const Teuchos::RCP<FieldEvaluator<GlobalOrdinal,SourceField> >& source_evaluator,
	Teuchos::RCP<FieldManager<TargetField> >& target_space_manager );

    //@{
    // Get the local indices of the target points that were not mapped.
    Teuchos::ArrayView<GlobalOrdinal>       getMissedTargetPoints();
    Teuchos::ArrayView<const GlobalOrdinal> getMissedTargetPoints() const;
    //@}

  private:

    // Compute globally unique ordinals for the target points.
    void computePointOrdinals( 
	const RCP_CoordFieldManager& target_coord_manager,
	Teuchos::Array<GlobalOrdinal>& target_ordinals );

    // Get the target points that are in the rendezvous decomposition box.
    void getTargetPointsInBox( 
	const BoundingBox& box,
	const CoordinateField& target_coords,
	const Teuchos::Array<GlobalOrdinal>& target_ordinals,
	Teuchos::Array<GlobalOrdinal>& targets_in_box );

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Map dimension.
    int d_dimension;

    // Boolean for storing missed points in the mapping.
    bool d_store_missed_points;

    // Process indexer for the source application.
    CommIndexer d_source_indexer;

    // Process indexer for the target application.
    CommIndexer d_target_indexer;

    // Indices for target points missed in the mapping.
    Teuchos::Array<GlobalOrdinal> d_missed_points;

    // Global-to-local ordinal map for target ordinals.
    std::tr1::unordered_map<GlobalOrdinal,GlobalOrdinal> d_target_g2l;

    // Source field map.
    RCP_TpetraMap d_source_map;

    // Target field map.
    RCP_TpetraMap d_target_map;

    // Source-to-target exporter.
    RCP_TpetraExport d_source_to_target_exporter;

    // Local source elements.
    Teuchos::Array<GlobalOrdinal> d_source_elements;

    // Local target coords.
    Teuchos::Array<double> d_target_coords;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_SharedDomainMap_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_SHAREDDOMAINMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_SharedDomainMap.hpp
//---------------------------------------------------------------------------//

