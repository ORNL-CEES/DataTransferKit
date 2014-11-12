//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_STKMeshEntityImpl_impl.hpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHENTITYIMPL_IMPL_HPP
#define DTK_STKMESHENTITYIMPL_IMPL_HPP

#include "DTK_DBC.hpp"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Bounding box extraction. Cartesian1D specialization.
template<>
void STKMeshEntityImpl::getNodeBounds( 
    const Teuchos::Array<stk::mesh::Entity>& entity_nodes,
    Teuchos::Tuple<double,6>& bounds,
    const Cartesian1DTag tag )
{
    int space_dim = 1;
    DTK_REQUIRE( physicalDimension() == space_dim );

    const stk::mesh::FieldBase* coord_field_base = 
	d_bulk_data->mesh_meta_data().coordinate_field();
    const stk::mesh::Field<double,stk::mesh::Cartesian1d>* coord_field =
	dynamic_cast<const stk::mesh::Field<double,stk::mesh::Cartesian1d>* >(
	    coord_field_base);

    Teuchos::Array<stk::mesh::Entity>::const_iterator entity_node_it;
    for ( entity_node_it = entity_nodes.begin();
	  entity_node_it != entity_nodes.end();
	  ++entity_node_it )
    {
	auto node_coords = stk::mesh::field_data( 
	    *coord_field, entity_node_it );
	for ( int d = 0; d < space_dim; ++d )
	{
	    bounds[d] = std::min( bounds[d], node_coords[d] );
	    bounds[d+3] = std::max( bounds[d+3], node_coords[d] );
	}
    }
}

//---------------------------------------------------------------------------//
// Bounding box extraction. Cartesian2D specialization.
template<>
void STKMeshEntityImpl::getNodeBounds( 
    const Teuchos::Array<stk::mesh::Entity>& entity_nodes,
    Teuchos::Tuple<double,6>& bounds,
    const Cartesian2DTag tag )
{
    int space_dim = 2;
    DTK_REQUIRE( physicalDimension() == space_dim );

    const stk::mesh::FieldBase* coord_field_base= 
	d_bulk_data->mesh_meta_data().coordinate_field();
    const stk::mesh::Field<double,stk::mesh::Cartesian1d>* coord_field =
	dynamic_cast<const stk::mesh::Field<double,stk::mesh::Cartesian2d>* >(
	    coord_field_base);

    Teuchos::Array<stk::mesh::Entity>::const_iterator entity_node_it;
    for ( entity_node_it = entity_nodes.begin();
	  entity_node_it != entity_nodes.end();
	  ++entity_node_it )
    {
	auto node_coords = stk::mesh::field_data( 
	    *coord_field, entity_node_it );
	for ( int d = 0; d < space_dim; ++d )
	{
	    bounds[d] = std::min( bounds[d], node_coords[d] );
	    bounds[d+3] = std::max( bounds[d+3], node_coords[d] );
	}
    }
}

//---------------------------------------------------------------------------//
// Bounding box extraction. Cartesian3D specialization.
template<>
void STKMeshEntityImpl::getNodeBounds( 
    const Teuchos::Array<stk::mesh::Entity>& entity_nodes,
    Teuchos::Tuple<double,6>& bounds,
    const Cartesian3DTag tag )
{
    int space_dim = 3;
    DTK_REQUIRE( physicalDimension() == space_dim );

    const stk::mesh::FieldBase* coord_field_base = 
	d_bulk_data->mesh_meta_data().coordinate_field();
    const stk::mesh::Field<double,stk::mesh::Cartesian1d>* coord_field =
	dynamic_cast<const stk::mesh::Field<double,stk::mesh::Cartesian3d>* >(
	    coord_field_base);

    Teuchos::Array<stk::mesh::Entity>::const_iterator entity_node_it;
    for ( entity_node_it = entity_nodes.begin();
	  entity_node_it != entity_nodes.end();
	  ++entity_node_it )
    {
	auto node_coords = stk::mesh::field_data( 
	    *coord_field, entity_node_it );
	for ( int d = 0; d < space_dim; ++d )
	{
	    bounds[d] = std::min( bounds[d], node_coords[d] );
	    bounds[d+3] = std::max( bounds[d+3], node_coords[d] );
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHENTITYIMPL_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityImpl_impl.hpp
//---------------------------------------------------------------------------//
