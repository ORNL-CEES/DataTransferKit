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
 * \brief Integral assembly map definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTEGRALASSEMBLYMAP_DEF_HPP
#define DTK_INTEGRALASSEMBLYMAP_DEF_HPP

#include <algorithm>
#include <limits>

#include "DTK_FieldTools.hpp"
#include "DTK_Assertion.hpp"
#include "DTK_Rendezvous.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Ptr.hpp>

#include <Tpetra_Import.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_MultiVector_def.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The communicator over which the map is generated.
 *
 * \param dimension The dimension of the map. This should be consistent with
 * all source and target objects (i.e. only 3 dimensional geometries will be
 * accepted with a 3 dimensional map). We need this here so we have a global
 * baseline for all objects that may or may not exist on all processes.
 */
template<class Mesh, class Geometry>
IntegralAssemblyMap<Mesh,Geometry>::IntegralAssemblyMap(
    const RCP_Comm& comm, const int dimension )
    : d_comm( comm )
    , d_dimension( dimension )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh, class Geometry>
IntegralAssemblyMap<Mesh,Geometry>::~IntegralAssemblyMap()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the integral map.
 *
 * \param source_mesh_manager Source mesh in the integral assembly problem. A
 * null RCP is a valid argument. This will be the case when a mesh manager is
 * only constructed on a subset of the processes that the integral assembly
 * map is constructed over. Note that the source mesh must exist only on
 * processes that reside within the IntegralAssemblyMap communicator.
 *
 * \param source_mesh_measure Measure interface object for getting mesh
 * element measures. A null RCP is a valid argument. This will be the case
 * when a mesh is only constructed on a subset of the processes that the
 * integral assembly map is constructed over.
 *
 * \param target_geometry_manager Target geometry in the integral assembly
 * problem. A null RCP is a valid argument. This will be the case when a
 * geometry manager is only constructed on a subset of the processes that the
 * integral assembly map is constructed over. Note that the target geometry
 * must exist only on processes that reside with the IntegralAssemblyMap
 * communicator.
 */
template<class Mesh, class Geometry>
void SharedDomainMap<Mesh,Geometry>::setup( 
    const RCP_MeshManager& source_mesh_manager, 
    const RCP_ElementMeasure& source_mesh_measure,
    const RCP_GeometryManager& target_geometry_manager )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the integral assembly map for a valid source field integrator and
 * target data space to the target points that were mapped.
 *
 * \param source_integrator Function integrator used to apply the mapping. This
 * FieldIntegrator must be valid for the source mesh used to generate the map.
 *
 * \param target_space_manager Target space into which the function
 * evaluations will be written. Enough space must be allocated to hold
 * evaluations at all points in all dimensions of the field.
 */
template<class Mesh, class Geometry>
template<class SourceField, class TargetField>
void IntegralAssemblyMap<Mesh,Geometry>::apply( 
    const Teuchos::RCP<FieldIntegrator<Mesh,SourceField> >& source_integrator,
    Teuchos::RCP<FieldManager<TargetField> >& target_space_manager )
{
    typedef FieldTraits<SourceField> SFT;
    typedef FieldTraits<TargetField> TFT;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute globally unique ordinals for the target points. Here an
 * invalid ordinal will be designated as the maximum value as specified by the
 * limits header for the ordinal type. We do this so that 0 may be a valid
 * ordinal.
 *
 * \param target_geometry_manager The geometries to compute global ordinals for.
 *
 * \param target_ordinals The computed globally unique ordinals for the target
 * geometries. 
 */
template<class Mesh, class Geometry>
void IntegralAssemblyMap<Mesh,Geometry>::computeGeometryOrdinals(
    const RCP_GeometryManager& target_geometry_manager,
    Teuchos::Array<GlobalOrdinal>& target_ordinals )
{
    // Set an existence value for the target geometry.
    bool target_exists = true;
    if ( target_geometry_manager.is_null() ) target_exists = false;
    int comm_rank = d_comm->getRank();
    GlobalOrdinal local_size = 0;

    if ( target_exists )
    {
	local_size = target_geometry_manager->numLocalGeometry();
    }
    d_comm->barrier();

    GlobalOrdinal global_max;
    Teuchos::reduceAll<int,GlobalOrdinal>( *d_comm,
					   Teuchos::REDUCE_MAX,
					   1,
					   &local_size,
					   &global_max );

    target_ordinals.resize( local_size );
    for ( GlobalOrdinal n = 0; n < local_size; ++n )
    {
	target_ordinals[n] = comm_rank*global_max + n;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_INTEGRALASSEMBLYMAP_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_IntegralAssemblyMap_def.hpp
//---------------------------------------------------------------------------//
