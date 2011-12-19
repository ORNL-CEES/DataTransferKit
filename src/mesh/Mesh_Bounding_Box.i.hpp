//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Mesh_Bounding_Box.i.hpp
 * \author Stuart Slattery
 * \date   Thu Dec 08 09:58:46 2011
 * \brief  Member definitions of class Mesh_Bounding_Box.
 */
//---------------------------------------------------------------------------//

#ifndef mesh_Mesh_Bounding_Box_i_hpp
#define mesh_Mesh_Bounding_Box_i_hpp

namespace Coupler
{

template<class HandleType, class CoordinateType>
bool Bounding_Box<HandleType,CoordinateType>::point_query(PointType point)
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

#endif // mesh_Mesh_Bounding_Box_i_hpp

//---------------------------------------------------------------------------//
//              end of mesh/Mesh_Bounding_Box.i.hpp
//---------------------------------------------------------------------------//
