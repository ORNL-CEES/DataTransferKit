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
 * \file DTK_MeshTools_def.hpp
 * \author Stuart R. Slattery
 * \brief MeshTools definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTOOLS_DEF_HPP
#define DTK_MESHTOOLS_DEF_HPP

#include <algorithm>
#include <iterator>

#include "DTK_DBC.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get a const view of the of the mesh vertex global ordinals. The
 * ArrayRCP object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A const view of the mesh vertex global ordinals. This view will not
 * manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<const typename MeshTools<Mesh>::GlobalOrdinal> 
MeshTools<Mesh>::verticesView( const Mesh& mesh )
{
    GlobalOrdinal num_vertices = std::distance( MT::verticesBegin( mesh ),
						MT::verticesEnd( mesh ) );

    if ( num_vertices == 0 )
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<const GlobalOrdinal>(
	    &*MT::verticesBegin(mesh), 0, num_vertices, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a non-const view of the of the mesh vertex global ordinals. The
 * ArrayRCP object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A non-const view of the mesh vertex global ordinals. The view will
 * not manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<typename MeshTools<Mesh>::GlobalOrdinal> 
MeshTools<Mesh>::verticesNonConstView( const Mesh& mesh )
{
    GlobalOrdinal num_vertices = std::distance( MT::verticesBegin( mesh ),
						MT::verticesEnd( mesh ) );

    if ( num_vertices == 0 )
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(
	    const_cast<GlobalOrdinal*>(&*MT::verticesBegin(mesh)), 0, 
	    num_vertices, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a view of the of the mesh vertex coordinates. The ArrayRCP
 * object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A const view of the mesh vertex coordinates. The view will not
 * manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<const double>
MeshTools<Mesh>::coordsView( const Mesh& mesh )
{
    GlobalOrdinal num_coords = std::distance( MT::coordsBegin( mesh ),
					      MT::coordsEnd( mesh ) );

    if ( num_coords == 0 )
    {
	return Teuchos::ArrayRCP<double>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<const double>(
	    &*MT::coordsBegin(mesh), 0, num_coords, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a non-const view of the of the mesh vertex coordinates. The
 * ArrayRCP object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A non-const view of the mesh vertex coordinates. The view will not
 * manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<double>
MeshTools<Mesh>::coordsNonConstView( const Mesh& mesh )
{
    GlobalOrdinal num_coords = std::distance( MT::coordsBegin( mesh ), 
					      MT::coordsEnd( mesh ) );

    if ( num_coords == 0 )
    {
	return Teuchos::ArrayRCP<double>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<double>(
	    const_cast<double*>(&*MT::coordsBegin(mesh)), 0, 
	    num_coords, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a view of the of the mesh element global ordinals. The ArrayRCP
 * object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A const view of the mesh element global ordinals. The view will not
 * manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<const typename MeshTools<Mesh>::GlobalOrdinal> 
MeshTools<Mesh>::elementsView( const Mesh& mesh )
{
    GlobalOrdinal num_elements = std::distance( MT::elementsBegin( mesh ), 
						MT::elementsEnd( mesh ) );
    if ( num_elements == 0 )
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<const GlobalOrdinal>(
	    &*MT::elementsBegin(mesh), 0, num_elements, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a non-const view of the of the mesh element global ordinals. The
 * ArrayRCP object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A non-const view of the mesh element global ordinals. The view will
 * not manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<typename MeshTools<Mesh>::GlobalOrdinal> 
MeshTools<Mesh>::elementsNonConstView( const Mesh& mesh )
{
    GlobalOrdinal num_elements = std::distance( MT::elementsBegin( mesh ), 
						MT::elementsEnd( mesh ) );
    
    if ( num_elements == 0 )
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(
	    const_cast<GlobalOrdinal*>(&*MT::elementsBegin(mesh)), 
	    0, num_elements, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a view of the of the mesh element connectivity. The ArrayRCP
 * object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A const view of the mesh element connectivity. The view will not
 * manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<const typename MeshTools<Mesh>::GlobalOrdinal> 
MeshTools<Mesh>::connectivityView( const Mesh& mesh )
{
    GlobalOrdinal num_connectivity = 
	std::distance( MT::connectivityBegin( mesh ),
		       MT::connectivityEnd( mesh ) );
    
    if ( num_connectivity == 0 )
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<const GlobalOrdinal>(
	    &*MT::connectivityBegin(mesh), 0, num_connectivity, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a non-const view of the of the mesh element connectivity. The
 * ArrayRCP object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A non-const view of the mesh element connectivity. The view will
 * not manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<typename MeshTools<Mesh>::GlobalOrdinal> 
MeshTools<Mesh>::connectivityNonConstView( const Mesh& mesh )
{
    GlobalOrdinal num_connectivity = 
	std::distance( MT::connectivityBegin( mesh ),
		       MT::connectivityEnd( mesh ) );
    
    if ( num_connectivity == 0 )
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<GlobalOrdinal>(
	    const_cast<GlobalOrdinal*>(&*MT::connectivityBegin(mesh)), 0, 
	    num_connectivity, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a view of the of the mesh connectivity permutation list. The
 * ArrayRCP object will not manage the memory. 
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A non-const view of the mesh connectivity permutation. The view
 * will not manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<const int> 
MeshTools<Mesh>::permutationView( const Mesh& mesh )
{
    int num_permutation = 
	std::distance( MT::permutationBegin( mesh ),
		       MT::permutationEnd( mesh ) );

    if ( num_permutation == 0 )
    {
	return Teuchos::ArrayRCP<int>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<const int>(
	    &*MT::permutationBegin( mesh ), 0, num_permutation, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a non-const view of the of the mesh connectivity permutation
 * list. The ArrayRCP object will not manage the memory.
 *
 * \param mesh The mesh to get the view for. This object must have mesh
 * traits. 
 *
 * \return A non-const view of the mesh connectivity permutation. The view
 * will not manage the memory.
 */
template <class Mesh> 
Teuchos::ArrayRCP<int> 
MeshTools<Mesh>::permutationNonConstView( const Mesh& mesh )
{
    int num_permutation = 
	std::distance( MT::permutationBegin( mesh ),
		       MT::permutationEnd( mesh ) );

    if ( num_permutation == 0 )
    {
	return Teuchos::ArrayRCP<int>(0,0);
    }
    else
    {
	return Teuchos::ArrayRCP<int>(
	    const_cast<int*>(&*MT::permutationBegin( mesh )), 0, 
	    num_permutation, false );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box for a mesh block.
 * 
 * \param mesh The mesh block for which to get the local bounding box.
 *
 * \return The local bounding box around the mesh block.
 */
template<class Mesh>
BoundingBox MeshTools<Mesh>::localBoundingBox( const Mesh& mesh )
{
    double x_min = -Teuchos::ScalarTraits<double>::rmax();
    double y_min = -Teuchos::ScalarTraits<double>::rmax();
    double z_min = -Teuchos::ScalarTraits<double>::rmax();

    double x_max = Teuchos::ScalarTraits<double>::rmax();
    double y_max = Teuchos::ScalarTraits<double>::rmax();
    double z_max = Teuchos::ScalarTraits<double>::rmax();

    GlobalOrdinal num_vertices = std::distance( MT::verticesBegin( mesh ),
						MT::verticesEnd( mesh ) );
    int vertex_dim = MT::vertexDim( mesh );

    if ( vertex_dim > 0 )
    {
	x_min = *std::min_element( 
	    MT::coordsBegin( mesh ),
	    MT::coordsBegin( mesh ) + num_vertices );
	x_max = *std::max_element( 
	    MT::coordsBegin( mesh ),
	    MT::coordsBegin( mesh ) + num_vertices );
    }
    if ( vertex_dim > 1 )
    {
	y_min = *std::min_element( 
	    MT::coordsBegin( mesh ) + num_vertices,
	    MT::coordsBegin( mesh ) + 2*num_vertices );
	y_max = *std::max_element( 
	    MT::coordsBegin( mesh ) + num_vertices,
	    MT::coordsBegin( mesh ) + 2*num_vertices );
    }
    if ( vertex_dim > 2 )
    {
	z_min = *std::min_element( 
	    MT::coordsBegin( mesh ) + 2*num_vertices,
	    MT::coordsBegin( mesh ) + 3*num_vertices );
	z_max = *std::max_element( 
	    MT::coordsBegin( mesh ) + 2*num_vertices,
	    MT::coordsBegin( mesh ) + 3*num_vertices );
    }

    return BoundingBox( x_min, y_min, z_min, x_max, y_max, z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global bounding box for a mesh block over the given
 * communicator.
 *
 * \param mesh The mesh block over which to get the global bounding box.
 *
 * \param comm The communicator over which the mesh block is defined.
 *
 * \return The global bounding box over the mesh block.
 */
template<class Mesh>
BoundingBox MeshTools<Mesh>::globalBoundingBox( const Mesh& mesh, 
						const RCP_Comm& comm )
{
    BoundingBox local_box = localBoundingBox( mesh );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();

    double global_x_min, global_y_min, global_z_min;
    double global_x_max, global_y_max, global_z_max;

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[0],
				    &global_x_min );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[1],
				    &global_y_min );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MIN,
				    1,
				    &local_bounds[2],
				    &global_z_min );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[3],
				    &global_x_max );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[4],
				    &global_y_max );

    Teuchos::reduceAll<int,double>( *comm, 
				    Teuchos::REDUCE_MAX,
				    1,
				    &local_bounds[5],
				    &global_z_max );

    return BoundingBox( global_x_min, global_y_min, global_z_min,
			global_x_max, global_y_max, global_z_max );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESHTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshTools_def.hpp
//---------------------------------------------------------------------------//

