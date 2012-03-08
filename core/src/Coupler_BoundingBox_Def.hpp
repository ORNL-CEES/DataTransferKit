//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Coupler_BoundingBox_Def.hpp
 * \author Stuart Slattery
 * \date   Thu Dec 08 09:58:46 2011
 * \brief  Member definitions of class BoundingBox.
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_BOUNDINGBOX_DEF_HPP
#define COUPLER_BOUNDINGBOX_DEF_HPP

namespace Coupler
{

template<class HandleType, class CoordinateType>
bool BoundingBox<HandleType,CoordinateType>::point_query(PointType point)
{
    bool return_val = false;

    if ( point.x() >= d_domain[0] &&
	 point.x() <= d_domain[1] &&
	 point.y() >= d_domain[2] &&
	 point.y() <= d_domain[3] &&
	 point.z() >= d_domain[4] &&
	 point.z() <= d_domain[5] )
    {
	return_val = true;
    }

    return return_val;
}

} // end namespace Coupler

#endif // COUPLER_BOUNDINGBOX_DEF_HPP

//---------------------------------------------------------------------------//
//              end of Coupler_BoundingBox_Def.hpp
//---------------------------------------------------------------------------//
