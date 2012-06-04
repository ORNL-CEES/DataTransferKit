//---------------------------------------------------------------------------//
/*!
 * \file DTK_BoundingBox.cpp
 * \author Stuart R. Slattery
 * \brief Bounding box definition.
 */
//---------------------------------------------------------------------------//

#include "DTK_BoundingBox.hpp"
#include <DTK_Exception.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
BoundingBox::BoundingBox()
    : d_x_min( 0.0 )
    , d_y_min( 0.0 )
    , d_z_min( 0.0 )
    , d_x_max( 0.0 )
    , d_y_max( 0.0 )
    , d_z_max( 0.0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
BoundingBox::BoundingBox( 
    const double x_min, const double y_min, const double z_min,
    const double x_max, const double y_max, const double z_max )
    : d_x_min( x_min )
    , d_y_min( y_min )
    , d_z_min( z_min )
    , d_x_max( x_max )
    , d_y_max( y_max )
    , d_z_max( z_max )
{
    testPrecondition( d_x_min < d_x_max, "Bounding box x_min > x_max" );
    testPrecondition( d_y_min < d_y_max, "Bounding box y_min > y_max" );
    testPrecondition( d_z_min < d_z_max, "Bounding box z_min > z_max" );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
BoundingBox::~BoundingBox()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \Determine if a point is in the box.
 */
bool BoundingBox::pointInBox( double coords[3] ) const
{
    if ( coords[0] >= d_x_min &&
	 coords[1] >= d_y_min &&
	 coords[2] >= d_z_min &&
	 coords[0] <= d_x_max &&
	 coords[1] <= d_y_max &&
	 coords[2] <= d_z_max )
    {
	return true;
    }

    return false;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BoundingBox.cpp
//---------------------------------------------------------------------------//

