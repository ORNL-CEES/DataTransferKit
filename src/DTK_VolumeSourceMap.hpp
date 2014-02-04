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
 * \file DTK_VolumeSourceMap.hpp
 * \author Stuart R. Slattery
 * \brief Volume source map declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_VOLUMESOURCEMAP_HPP
#define DTK_VOLUMESOURCEMAP_HPP

#include <boost/tr1/unordered_map.hpp>

#include "DTK_GeometryTraits.hpp"
#include "DTK_GeometryManager.hpp"
#include "DTK_FieldTraits.hpp"
#include "DTK_FieldEvaluator.hpp"
#include "DTK_FieldManager.hpp"
#include "DTK_CommIndexer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Directory.hpp>
#include <Tpetra_Import.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class VolumeSourceMap
 * \brief A map for transfers where the volume is a source.
 *
 * In all cases supported by this map, the source consists of a group of
 * volumes distributed in parallel. That target will provide a list of
 * coordinates that data is desired for. These may be the centroids of target
 * volumes, quadrature points in a target mesh, or any other representation
 * that may be resolved by points. The apply stage will provide the source
 * with a list of local source geometries and a set of target points at which
 * to evaluate the source function. This evaluation can be interpreted in any
 * way, simply moving the volume data, doing a zero order evaluation for the
 * volume to volume case, or a higher order functional evaluation for the
 * quadrature point case.
 */
//---------------------------------------------------------------------------//
template<class Geometry, class GlobalOrdinal, class CoordinateField>
class VolumeSourceMap
{
  public:

    //@{
    //! Typedefs.
    typedef Geometry                                  geometry_type;
    typedef GeometryTraits<Geometry>                  GT;
    typedef GeometryManager<Geometry,GlobalOrdinal>   GeometryManagerType;
    typedef Teuchos::RCP<GeometryManagerType>         RCP_GeometryManager;
    typedef GlobalOrdinal                             global_ordinal_type;
    typedef CoordinateField                           coord_field_type;
    typedef FieldTraits<CoordinateField>              CFT;
    typedef typename CFT::size_type                   CoordOrdinal;
    typedef FieldManager<CoordinateField>             CoordFieldManagerType;
    typedef Teuchos::RCP<CoordFieldManagerType>       RCP_CoordFieldManager;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;
    typedef Tpetra::Map<int,GlobalOrdinal>            TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    typedef Tpetra::Import<int,GlobalOrdinal>         ImportType;
    typedef Teuchos::RCP<ImportType>                  RCP_TpetraImport;
    //@}

    // Constructor.
    VolumeSourceMap( const RCP_Comm& comm, const int dimension,
		     bool store_missed_points = false,
		     const double geometric_tolerance = 1.0e-6 );

    // Destructor.
    ~VolumeSourceMap();

    // Generate the volume source map.
    void setup( const RCP_GeometryManager& source_geometry_manager, 
		const RCP_CoordFieldManager& target_coord_manager );

    //@{
    // Get the local indices of the target points that were not mapped.
    Teuchos::ArrayView<GlobalOrdinal>       getMissedTargetPoints();
    Teuchos::ArrayView<const GlobalOrdinal> getMissedTargetPoints() const;
    //@}

    // Apply the volume source map by evaluating a function at the target points
    // that were mapped.
    template<class SourceField, class TargetField>
    void apply( 
	const Teuchos::RCP<FieldEvaluator<GlobalOrdinal,SourceField> >& source_evaluator,
	Teuchos::RCP<FieldManager<TargetField> >& target_space_manager );

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

    // Geometric tolerance.
    double d_geometric_tolerance;

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

    // Source-to-target importer.
    RCP_TpetraImport d_source_to_target_importer;

    // Local source geometries.
    Teuchos::Array<GlobalOrdinal> d_source_geometry;

    // Local target coords.
    Teuchos::Array<double> d_target_coords;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_VolumeSourceMap_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_VOLUMESOURCEMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_VolumeSourceMap.hpp
//---------------------------------------------------------------------------//

