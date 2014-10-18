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

#include <boost/tr1/unordered_map.hpp>

#include "DTK_MeshTypes.hpp"
#include "DTK_MeshManager.hpp"
#include "DTK_GeometryTraits.hpp"
#include "DTK_GeometryManager.hpp"
#include "DTK_ElementMeasure.hpp"
#include "DTK_FieldIntegrator.hpp"
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
  \class IntegralAssemblyMap
  \brief A map for assembling integrals over geometric objects using a
  conformal mesh assumption.

  Here we assume that the mesh supplied to the map is conformal to the
  geometric objects provided by the target. For example, if the target
  geometry consists of a group of cylinders, the source mesh should include
  the mesh conformal to those cylinders. 
*/
template<class Geometry>
class IntegralAssemblyMap
{
  public:

    //@{
    //! Typedefs.
    typedef Geometry                                  geometry_type;
    typedef GeometryTraits<Geometry>                  GT;
    typedef GeometryManager<Geometry,MeshId>          GeometryManagerType;
    typedef Teuchos::RCP<GeometryManagerType>         RCP_GeometryManager;
    //@}

    // Constructor.
    IntegralAssemblyMap( const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
			 const int dimension, 
			 const double geometric_tolerance = 1.0e-6, 
			 bool all_vertices_for_inclusion = true );

    // Destructor.
    ~IntegralAssemblyMap();

    // Generate the integral assembly map.
    void setup( 
	const Teuchos::RCP<MeshManager>& source_mesh_manager,
	const Teuchos::RCP<ElementMeasure>& source_mesh_measure,
	const RCP_GeometryManager& target_geometry_manager );

    // Apply the shared domain map by integrating the source function over the
    // target geometries.
    template<class SourceField, class TargetField>
    void apply(
	const Teuchos::RCP<FieldIntegrator<SourceField> >& source_integrator,
	Teuchos::RCP<FieldManager<TargetField> >& target_space_manager );

  private:

    // Compute globally unique ordinals for the target geometries.
    void computeGeometryOrdinals( 
	const RCP_GeometryManager& target_geometry_manager,
	Teuchos::Array<MeshId>& target_ordinals );

  private:

    // Communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Map dimension.
    int d_dimension;

    // Geometric tolerance.
    double d_geometric_tolerance;
    
    // Flag for element-in-geometry vertex inclusion requirement.
    bool d_all_vertices_for_inclusion;

    // Process indexer for the source application.
    CommIndexer d_source_indexer;

    // Process indexer for the target application.
    CommIndexer d_target_indexer;

    // Source field map.
    Teuchos::RCP<const Tpetra::Map<int,MeshId> > d_source_map;

    // Target field map.
    Teuchos::RCP<const Tpetra::Map<int,MeshId> > d_target_map;

    // Source-to-target importer.
    Teuchos::RCP<const Tpetra::Import<int,MeshId> > d_source_to_target_importer;

    // Local source elements to drive the function integrations (source
    // decomposition).
    Teuchos::Array<MeshId> d_source_elements;

    // Local geometry measures from element measures.
    Teuchos::Array<double> d_geometry_measures;

    // Source element local ordinals that construct the local geometry
    // integral. 
    Teuchos::Array<Teuchos::Array<MeshId> > d_integral_elements;

    // Global-to-local ordinal map for the target geometry in the target
    // decomposition. 
    std::tr1::unordered_map<
	MeshId, typename Teuchos::ArrayRCP<Geometry>::size_type> d_target_g2l;
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

