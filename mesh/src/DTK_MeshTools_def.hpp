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

#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box for a mesh.
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

    typename MT::global_ordinal_type num_nodes =
	std::distance( MT::nodesBegin( mesh ),
		       MT::nodesEnd( mesh ) );
    std::size_t node_dim = MT::nodeDim( mesh );

    if ( node_dim > 0 )
    {
	x_min = *std::min_element( 
	    MT::coordsBegin( mesh ),
	    MT::coordsBegin( mesh ) + num_nodes );
	x_max = *std::max_element( 
	    MT::coordsBegin( mesh ),
	    MT::coordsBegin( mesh ) + num_nodes );
    }
    if ( node_dim > 1 )
    {
	y_min = *std::min_element( 
	    MT::coordsBegin( mesh ) + num_nodes,
	    MT::coordsBegin( mesh ) + 2*num_nodes );
	y_max = *std::max_element( 
	    MT::coordsBegin( mesh ) + num_nodes,
	    MT::coordsBegin( mesh ) + 2*num_nodes );
    }
    if ( node_dim > 2 )
    {
	z_min = *std::min_element( 
	    MT::coordsBegin( mesh ) + 2*num_nodes,
	    MT::coordsBegin( mesh ) + 3*num_nodes );
	z_max = *std::max_element( 
	    MT::coordsBegin( mesh ) + 2*num_nodes,
	    MT::coordsBegin( mesh ) + 3*num_nodes );
    }

    return BoundingBox( x_min, y_min, z_min, x_max, y_max, z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global bounding box for a mesh.
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

