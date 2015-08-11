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
 * \brief DTK_ClassicMeshElementLocalMap_impl.hpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSICMESHELEMENTLOCALMAP_IMPL_HPP
#define DTK_CLASSICMESHELEMENTLOCALMAP_IMPL_HPP

#include "DTK_ClassicMeshElementExtraData.hpp"
#include "DTK_IntrepidCellLocalMap.hpp"
#include "DTK_DBC.hpp"

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Mesh>
ClassicMeshElementLocalMap<Mesh>::ClassicMeshElementLocalMap(
    const Teuchos::RCP<ClassicMesh<Mesh> >& mesh )
    : d_mesh( mesh )
    , d_inclusion_tol( 1.0e-6 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Set parameters for mapping.
template<class Mesh>
void ClassicMeshElementLocalMap<Mesh>::setParameters(
    const Teuchos::ParameterList& parameters )
{
    if ( parameters.isParameter("Point Inclusion Tolerance") )
    {	    
	d_inclusion_tol = parameters.get<double>("Point Inclusion Tolerance");
    }
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
template<class Mesh>
double ClassicMeshElementLocalMap<Mesh>::measure( const Entity& entity ) const
{
    // Get the block id and topology.
    int block_id = Teuchos::rcp_dynamic_cast<ClassicMeshElementExtraData>(
	entity.extraData() )->d_block_id;
    shards::CellTopology entity_topo = d_mesh->getBlockTopology( block_id );

    // Get the entity coordinates.
    Intrepid::FieldContainer<double> entity_coords =
	d_mesh->getElementNodeCoordinates( entity.id(), block_id );

    // Compute the measure.
    return IntrepidCellLocalMap::measure( entity_topo, entity_coords );
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
template<class Mesh>
void ClassicMeshElementLocalMap<Mesh>::centroid( 
    const Entity& entity, const Teuchos::ArrayView<double>& centroid ) const
{ 
    // Get the block id and topology.
    int block_id = Teuchos::rcp_dynamic_cast<ClassicMeshElementExtraData>(
	entity.extraData() )->d_block_id;
    shards::CellTopology entity_topo = d_mesh->getBlockTopology( block_id );

    // Get the entity coordinates.
    Intrepid::FieldContainer<double> entity_coords =
	d_mesh->getElementNodeCoordinates( entity.id(), block_id );

    // Compute the centroid of the element.
    IntrepidCellLocalMap::centroid( entity_topo, entity_coords, centroid );
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance.
template<class Mesh>
bool ClassicMeshElementLocalMap<Mesh>::isSafeToMapToReferenceFrame(
    const Entity& entity,
    const Teuchos::ArrayView<const double>& physical_point ) const
{
    return EntityLocalMap::isSafeToMapToReferenceFrame(
	entity, physical_point );
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
template<class Mesh>
bool ClassicMeshElementLocalMap<Mesh>::mapToReferenceFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& physical_point,
    const Teuchos::ArrayView<double>& reference_point ) const
{ 
    // Get the block id and topology.
    int block_id = Teuchos::rcp_dynamic_cast<ClassicMeshElementExtraData>(
	entity.extraData() )->d_block_id;
    shards::CellTopology entity_topo = d_mesh->getBlockTopology( block_id );

    // Get the entity coordinates.
    Intrepid::FieldContainer<double> entity_coords =
	d_mesh->getElementNodeCoordinates( entity.id(), block_id );

    // Use the cell to perform the element mapping.
    return IntrepidCellLocalMap::mapToReferenceFrame(
	entity_topo, entity_coords, physical_point, reference_point );
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
template<class Mesh>
bool ClassicMeshElementLocalMap<Mesh>::checkPointInclusion( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    // Get the block id and topology.
    int block_id = Teuchos::rcp_dynamic_cast<ClassicMeshElementExtraData>(
	entity.extraData() )->d_block_id;
    shards::CellTopology entity_topo = d_mesh->getBlockTopology( block_id );

    // Get the entity coordinates.
    Intrepid::FieldContainer<double> entity_coords =
	d_mesh->getElementNodeCoordinates( entity.id(), block_id );

    // Check point inclusion in the element.
    return IntrepidCellLocalMap::checkPointInclusion( 
	entity_topo, reference_point, d_inclusion_tol );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
template<class Mesh>
void ClassicMeshElementLocalMap<Mesh>::mapToPhysicalFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& physical_point ) const
{
    // Get the block id and topology.
    int block_id = Teuchos::rcp_dynamic_cast<ClassicMeshElementExtraData>(
	entity.extraData() )->d_block_id;
    shards::CellTopology entity_topo = d_mesh->getBlockTopology( block_id );

    // Get the entity coordinates.
    Intrepid::FieldContainer<double> entity_coords =
	d_mesh->getElementNodeCoordinates( entity.id(), block_id );

    // Map from the element.
    IntrepidCellLocalMap::mapToPhysicalFrame( 
	entity_topo, entity_coords, reference_point, physical_point );
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference point.
template<class Mesh>
void ClassicMeshElementLocalMap<Mesh>::normalAtReferencePoint( 
    const Entity& entity,
    const Entity& parent_entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& normal ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end  DTK_CLASSICMESHELEMENTLOCALMAP_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicMeshElementLocalMap_impl.hpp
//---------------------------------------------------------------------------//
