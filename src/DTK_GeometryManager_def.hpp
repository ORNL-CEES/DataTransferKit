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

#include "DTK_Assertion.hpp"
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
 * \param comm The communicator over which the geometry is defined.
 *
 * \param dim The dimension of the geometry.
 */
template<class Geometry>
GeometryManager<Geometry>::GeometryManager( 
    const Teuchos::ArrayRCP<Geometry>& geometry, const RCP_Comm& comm,
    const int dim )
    : d_geometry( geometry )
    , d_comm( comm )
    , d_dim( dim )
{
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
template<class Geometry>
GeometryManager<Geometry>::~GeometryManager()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global number of objects owned by this manager.
 *
 * \return The global number of objects owned by this manager.
 */
template<class Geometry>
const typename Teuchos::ArrayRCP<Geometry>::size_type
GeometryManager<Geometry>::globalNumGeometry() const
{
    typename Teuchos::ArrayRCP<Geometry>::size_type global_size =
	d_geometry.size();

    Teuchos::reduceAll<int,typename Teuchos::ArrayRCP<Geometry>::size_type>( 
	*d_comm, Teuchos::REDUCE_SUM, 1, &global_size, &global_size );

    return global_size;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the bounding boxes for the objects owned by this manager.
 *
 * \return The bounding boxes for the geometry owned by this manager.
 */
template<class Geometry>
Teuchos::Array<BoundingBox> GeometryManager<Geometry>::boundingBoxes() const
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
 * \brief Validate the geometry to the domain model.
 */
template<class Geometry>
void GeometryManager<Geometry>::validate()
{
    // Dimensions greater than 3 are not valid.
    testPrecondition( 0 <= d_dim && d_dim <= 3 );

    // Check that all local geometries have the same dimension.
    typename Teuchos::ArrayRCP<Geometry>::const_iterator geom_iterator;
    for ( geom_iterator = d_geometry.begin();
	  geom_iterator != d_geometry.end();
	  ++geom_iterator )
    {
	testPrecondition( GT::dim( *geom_iterator ) == d_dim );
    }

    // Check that the geometry dimension is the same on every node.
    Teuchos::Array<int> local_dims( d_comm->getSize(), 0 );
    local_dims[ d_comm->getRank() ] = d_dim;
    Teuchos::reduceAll<int,int>( *d_comm, Teuchos::REDUCE_SUM,
				 local_dims.size(),
				 &local_dims[0], &local_dims[0] ); 
    Teuchos::Array<int>::iterator unique_bound;
    unique_bound = std::unique( local_dims.begin(), local_dims.end() );
    int unique_dim = std::distance( local_dims.begin(), unique_bound );
    testPrecondition( 1 == unique_dim );
    local_dims.clear();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_GEOMETRYMANAGER_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryManager_def.hpp
//---------------------------------------------------------------------------//

