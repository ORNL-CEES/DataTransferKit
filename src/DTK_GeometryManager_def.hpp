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
 * \file DTK_GeometryManager_def.hpp
 * \author Stuart R. Slattery
 * \brief Geometry manager definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYMANAGER_DEF_HPP
#define DTK_GEOMETRYMANAGER_DEF_HPP

#include "DTK_DBC.hpp"
#include "DataTransferKit_config.hpp"

#include <Teuchos_CommHelpers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor. If Design-By-Contract is enabled, the constructor will
 * validate the geometry description to the domain model. This requires a few
 * global communications.
 *
 * \param geometry The geometry that this object is managing. This geometry
 * must have GeometryTraits.
 * 
 * \param geom_gids The global ordinals for the geometry this object is
 * managing. 
 *
 * \param comm The communicator over which the geometry is defined.
 *
 * \param dim The dimension of the geometry.
 */
template<class Geometry,class GlobalOrdinal>
GeometryManager<Geometry,GlobalOrdinal>::GeometryManager( 
    const Teuchos::ArrayRCP<Geometry>& geometry,
    const Teuchos::ArrayRCP<GlobalOrdinal>& geom_gids,
    const RCP_Comm& comm, const int dim )
    : d_geometry( geometry )
    , d_geom_gids( geom_gids )
    , d_comm( comm )
    , d_dim( dim )
{
    DTK_REQUIRE( d_geometry.size() == d_geom_gids.size() );

    // If we're checking with Design-by-Contract, validate the geometry to the
    // domain model.
#if HAVE_DTK_DBC
    validate();
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Geometry,class GlobalOrdinal>
GeometryManager<Geometry,GlobalOrdinal>::~GeometryManager()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global number of objects owned by this manager.
 *
 * \return The global number of objects owned by this manager.
 */
template<class Geometry,class GlobalOrdinal>
const typename Teuchos::ArrayRCP<Geometry>::size_type
GeometryManager<Geometry,GlobalOrdinal>::globalNumGeometry() const
{
    typename Teuchos::ArrayRCP<Geometry>::size_type global_size =
	d_geometry.size();
    typename Teuchos::ArrayRCP<Geometry>::size_type global_size_copy =
	d_geometry.size();

    Teuchos::reduceAll<int,typename Teuchos::ArrayRCP<Geometry>::size_type>( 
	*d_comm, Teuchos::REDUCE_SUM, 1, &global_size, &global_size_copy );

    return global_size_copy;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the bounding boxes for the objects owned by this manager.
 *
 * \return The bounding boxes for the geometry owned by this manager.
 */
template<class Geometry,class GlobalOrdinal>
Teuchos::Array<BoundingBox> 
GeometryManager<Geometry,GlobalOrdinal>::boundingBoxes() const
{
    Teuchos::Array<BoundingBox> boxes( d_geometry.size() );
    Teuchos::Array<BoundingBox>::iterator box_iterator;
    typename Teuchos::ArrayRCP<Geometry>::const_iterator geom_iterator;
    for ( geom_iterator = d_geometry.begin(), box_iterator = boxes.begin();
	  geom_iterator != d_geometry.end();
	  ++geom_iterator, ++box_iterator )
    {
	*box_iterator = GT::boundingBox( *geom_iterator );
    }

    return boxes;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box for the objects owned by this manager.
 */
template<class Geometry,class GlobalOrdinal>
BoundingBox GeometryManager<Geometry,GlobalOrdinal>::localBoundingBox() const
{
    double global_x_min = Teuchos::ScalarTraits<double>::rmax();
    double global_y_min = Teuchos::ScalarTraits<double>::rmax();
    double global_z_min = Teuchos::ScalarTraits<double>::rmax();
    double global_x_max = -Teuchos::ScalarTraits<double>::rmax();
    double global_y_max = -Teuchos::ScalarTraits<double>::rmax();
    double global_z_max = -Teuchos::ScalarTraits<double>::rmax();

    // Get the local bounding boxes compute the local bounding box.
    BoundingBox local_box;
    Teuchos::Tuple<double,6> box_bounds;
    Teuchos::Array<BoundingBox> boxes = boundingBoxes();
    Teuchos::Array<BoundingBox>::const_iterator box_iterator;
    DTK_CHECK( !boxes.empty() );
    for ( box_iterator = boxes.begin();
	  box_iterator != boxes.end();
	  ++box_iterator )
    {
	box_bounds = box_iterator->getBounds();

	if ( box_bounds[0] < global_x_min )
	{
	    global_x_min = box_bounds[0];
	}
	if ( box_bounds[1] < global_y_min )
	{
	    global_y_min = box_bounds[1];
	}
	if ( box_bounds[2] < global_z_min )
	{
	    global_z_min = box_bounds[2];
	}
	if ( box_bounds[3] > global_x_max )
	{
	    global_x_max = box_bounds[3];
	}
	if ( box_bounds[4] > global_y_max )
	{
	    global_y_max = box_bounds[4];
	}
	if ( box_bounds[5] > global_z_max )
	{
	    global_z_max = box_bounds[5];
	}
    }

    return BoundingBox( global_x_min, global_y_min, global_z_min,
			global_x_max, global_y_max, global_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global bounding box for the objects owned by this manager.
 */
template<class Geometry,class GlobalOrdinal>
BoundingBox GeometryManager<Geometry,GlobalOrdinal>::globalBoundingBox() const
{
    Teuchos::Tuple<double,6> local_bounds =
	Teuchos::tuple( Teuchos::ScalarTraits<double>::rmax(),
			Teuchos::ScalarTraits<double>::rmax(),
			Teuchos::ScalarTraits<double>::rmax(),
			-Teuchos::ScalarTraits<double>::rmax(),
			-Teuchos::ScalarTraits<double>::rmax(),
			-Teuchos::ScalarTraits<double>::rmax() );
    if ( d_geometry.size() > 0 )
    {
	BoundingBox local_box = localBoundingBox();
	local_bounds = local_box.getBounds();
    }

    double global_x_min, global_y_min, global_z_min;
    double global_x_max, global_y_max, global_z_max;

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[0],
				    &global_x_min );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[1],
				    &global_y_min );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[2],
				    &global_z_min );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[3],
				    &global_x_max );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[4],
				    &global_y_max );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[5],
				    &global_z_max );

    return BoundingBox( global_x_min, global_y_min, global_z_min,
			global_x_max, global_y_max, global_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Validate the geometry to the domain model.
 */
template<class Geometry,class GlobalOrdinal>
void GeometryManager<Geometry,GlobalOrdinal>::validate()
{
    // Dimensions greater than 3 are not valid.
    DTK_REQUIRE( 0 <= d_dim && d_dim <= 3 );

    // Check that all local geometries have the same dimension.
    typename Teuchos::ArrayRCP<Geometry>::const_iterator geom_iterator;
    for ( geom_iterator = d_geometry.begin();
	  geom_iterator != d_geometry.end();
	  ++geom_iterator )
    {
	DTK_REQUIRE( GT::dim( *geom_iterator ) == d_dim );
    }

    // Check that the geometry dimension is the same on every node.
    Teuchos::Array<int> local_dims( d_comm->getSize(), 0 );
    Teuchos::Array<int> local_dims_copy( d_comm->getSize(), 0 );
    local_dims[ d_comm->getRank() ] = d_dim;
    Teuchos::reduceAll<int,int>( *d_comm, Teuchos::REDUCE_SUM,
				 local_dims.size(),
				 &local_dims[0], &local_dims_copy[0] ); 
    Teuchos::Array<int>::iterator unique_bound;
    std::sort( local_dims_copy.begin(), local_dims_copy.end() );
    unique_bound = 
	std::unique( local_dims_copy.begin(), local_dims_copy.end() );
    int unique_dim = std::distance( local_dims_copy.begin(), unique_bound );
    DTK_REQUIRE( 1 == unique_dim );
    local_dims_copy.clear();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_GEOMETRYMANAGER_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryManager_def.hpp
//---------------------------------------------------------------------------//

