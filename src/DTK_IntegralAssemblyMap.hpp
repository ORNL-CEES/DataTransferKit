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
 * \file DTK_IntegralAssemblyMap.hpp
 * \author Stuart R. Slattery
 * \brief Integral assembly map declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTEGRALASSEMBLYMAP_HPP
#define DTK_INTEGRALASSEMBLYMAP_HPP

#include "DTK_MeshTraits.hpp"
#include "DTK_MeshManager.hpp"
#include "DTK_GeometryTraits.hpp"
#include "DTK_GeometryManager.hpp"
#include "DTK_ElementMeasure.hpp"
#include "DTK_FieldIntegrator.hpp"
#include "DTK_BoundingBox.hpp"
#include "DTK_CommIndexer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Map_def.hpp>
#include <Tpetra_Directory_decl.hpp>
#include <Tpetra_Directory_def.hpp>
#include <Tpetra_Export.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class IntegralAssemblyMap
  \brief A map for assembling integrals over geometric objects using a
  conformal mesh assumption.

  Here we assume that the mesh supplied to the map is conformal to the
  geometric objects provided by the target For example, if the target
  geometry consists of a group of cylinders, the source mesh should be the
  mesh over those cylinders and only that mesh, not the mesh in the domain
  between the cylinders.
*/
template<class Mesh, class Geometry>
class IntegralAssemblyMap
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
    typedef ElementMeasure<Mesh>                      ElementMeasureType;
    typedef Teuchos::RCP<ElementMeasureType>          RCP_ElementMeasure;
    typedef Geometry                                  geometry_type;
    typedef GeometryTraits<Geometry>                  GT;
    typedef GeometryManager<Geometry>                 GeometryManagerType;
    typedef Teuchos::RCP<GeometryManagerType>         RCP_GeometryManager;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;
    typedef Tpetra::Map<GlobalOrdinal>                TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    typedef Tpetra::Export<GlobalOrdinal>             ExportType;
    typedef Teuchos::RCP<ExportType>                  RCP_TpetraExport;
    //@}

    // Constructor.
    IntegralAssemblyMap( const RCP_Comm& comm, const int dimension );

    // Destructor.
    ~IntegralAssemblyMap();

    // Generate the integral assembly map.
    void setup( 
	const RCP_MeshManager& source_mesh_manager,
	const RCP_ElementMeasure& source_mesh_measure,
	const RCP_GeometryManager& target_geometry_manager );

    // Apply the shared domain map by integrating the source function over the
    // target geometries.
    template<class SourceField, class TargetField>
    void apply(
	const Teuchos::RCP<FieldIntegrator<Mesh,SourceField> >& source_integrator,
	Teuchos::RCP<FieldManager<TargetField> >& target_space_manager );

  private:

    // Compute globally unique ordinals for the target geometries.
    void computeGeometryOrdinals( 
	const RCP_GeometryManager& target_geometry_manager,
	Teuchos::Array<GlobalOrdinal>& target_ordinals );

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Map dimension.
    int d_dimension;

    // Process indexer for the source application.
    CommIndexer d_source_indexer;

    // Process indexer for the target application.
    CommIndexer d_target_indexer;

    // Global-to-local ordinal map for target ordinals.
    std::map<GlobalOrdinal,GlobalOrdinal> d_target_g2l;

    // Source field map.
    RCP_TpetraMap d_source_map;

    // Target field map.
    RCP_TpetraMap d_target_map;

    // Source-to-target exporter.
    RCP_TpetraExport d_source_to_target_exporter;

    // Local source elements (source decomposition).
    Teuchos::Array<GlobalOrdinal> d_source_elements;

    // Local source element measures (target decomposition).
    Teuchos::Array<double> d_source_element_measures;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_IntegralAssemblyMap_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_INTEGRALASSEMBLYMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_IntegralAssemblyMap.hpp
//---------------------------------------------------------------------------//

