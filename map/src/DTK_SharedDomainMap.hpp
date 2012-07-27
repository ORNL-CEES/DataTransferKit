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

#include <DTK_FieldTraits.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Export.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class SharedDomainMap
 * \brief A map for shared domain problems.
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
    typedef typename MeshManagerType::BlockIterator   BlockIterator;
    typedef CoordinateField                           coord_field_type;
    typedef FieldTraits<CoordinateField>              CFT;
    typedef typename CFT::size_type                   CoordOrdinal;
    typedef FieldManager<CoordinateField>             FieldManagerType;
    typedef Teuchos::RCP<FieldManagerType>            RCP_FieldManager;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;
    typedef Tpetra::Map<GlobalOrdinal>                TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    typedef Tpetra::Export<GlobalOrdinal>             ExportType;
    typedef Teuchos::RCP<ExportType>                  RCP_Export;
    //!@}

    // Constructor.
    SharedDomainMap( const RCP_Comm& comm, bool keep_missed_points = false );

    // Destructor.
    ~SharedDomainMap();

    // Generate the shared domain map.
    void setup( const RCP_MeshManager& source_mesh_manager, 
		const RCP_FieldManager& target_coord_manager );

    // Get the local indices of the target points that were not mapped.
    Teuchos::ArrayView<const CoordOrdinal> getMissedTargetPoints() const;

    // Apply the shared domain map to the target points that were mapped.
    template<class SourceField, class TargetField>
    void apply( const Teuchos::RCP< FieldEvaluator<Mesh,SourceField> >& 
		source_evaluator,
		Teuchos::RCP< FieldManager<TargetField> >& target_space_manager );

  private:

    // Get the target points that are in the rendezvous decomposition box.
    void getTargetPointsInBox( const BoundingBox& box,
			       const CoordinateField& target_coords,
			       Teuchos::Array<CoordOrdinal>& points_in_box );

    // Compute globally unique ordinals for the target points.
    void computePointOrdinals( const CoordinateField& target_coords,
			       Teuchos::Array<GlobalOrdinal>& ordinals );

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Boolean for keeping missed points in the mapping.
    bool d_keep_missed_points;

    // Indices for target points missed in the mapping.
    Teuchos::Array<CoordOrdinal> d_missed_points;

    // Export field map.
    RCP_TpetraMap d_export_map;

    // Import field map.
    RCP_TpetraMap d_import_map;

    // Field data exporter.
    RCP_Export d_data_export;

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

