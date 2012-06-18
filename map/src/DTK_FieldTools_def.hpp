//---------------------------------------------------------------------------//
/*!
 * \file DTK_FieldTools_def.hpp
 * \author Stuart R. Slattery
 * \brief FieldTools definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDTOOLS_DEF_HPP
#define DTK_FIELDTOOLS_DEF_HPP

#include <algorithm>

#include <DTK_Exception.hpp>

#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box for a coordinate field.
 */
template<class Field>
BoundingBox FieldTools<Field>::coordLocalBoundingBox( const Field& field )
{
    double x_min = Teuchos::ScalarTraits<double>::rmin();
    double y_min = Teuchos::ScalarTraits<double>::rmin();
    double z_min = Teuchos::ScalarTraits<double>::rmin();

    double x_max = Teuchos::ScalarTraits<double>::rmax();
    double y_max = Teuchos::ScalarTraits<double>::rmax();
    double z_max = Teuchos::ScalarTraits<double>::rmax();

    std::size_t dim = FT::dim( field );

    FT::size_type num_elements = std::distance( FT::begin( field ),
						FT::end( field ) ) / dim;


    if ( dim > 0 )
    {
	x_min = *std::min_element( 
	    FT::begin( field ),
	    FT::begin( field ) + num_elements );
	x_max = *std::max_element( 
	    FT::begin( field ),
	    FT::begin( field ) + num_elements );
    }
    if ( dim > 1 )
    {
	y_min = *std::min_element( 
	    FT::begin( field ) + num_elements,
	    FT::begin( field ) + 2*num_elements );
	y_max = *std::max_element( 
	    FT::begin( field ) + num_elements,
	    FT::begin( field ) + 2*num_elements );
    }
    if ( dim > 2 )
    {
	z_min = *std::min_element( 
	    FT::begin( field ) + 2*num_elements,
	    FT::begin( field ) + 3*num_elements );
	z_max = *std::max_element( 
	    FT::begin( field ) + 2*num_elements,
	    FT::begin( field ) + 3*num_elements );
    }
    if ( dim > 3 )
    {
	throw InvariantException( 
	    "Nodes with greater than 3 dimensions not supported" );
    }

    return BoundingBox( x_min, y_min, z_min, x_max, y_max, z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global bounding box for a coordinate field.
 */
template<class Field>
BoundingBox FieldTools<Field>::coordGlobalBoundingBox( const Field& field,
						       const RCP_Comm& comm )
{
    BoundingBox local_box = coordLocalBoundingBox( field );
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

#endif // end DTK_FIELDTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldTools_def.hpp
//---------------------------------------------------------------------------//

