//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Point.hpp
 * \author Stuart R. Slattery
 * \date   Wed May 25 12:23:57 2011
 * \brief  Cartesian Point class definition
 */
//---------------------------------------------------------------------------//

#ifndef mesh_Point_hpp
#define mesh_Point_hpp

namespace Coupler
{
//===========================================================================//
/*!
 * \class Point
 * \brief Cartesian Point class container. Holds (x,y,z)
 *
 * \par Code Sample:
 * \code
 *     double x = 0.0;
 *     double y = 1.0;
 *     double z = 1.0;
 *     Point pt(x,y,z);
 * \endcode
 */
/*! 
 * \example mesh_type/test/tstPoint.cc
 *
 * Test of Point.
 */
//===========================================================================//

template<class HandleType_T, class CoordinateType_T>
class Point 
{
  public:

    //@{ 
    //! Useful typedefs.
    typedef HandleType_T                  HandleType;
    typedef CoordinateType_T              CoordinateType;
    //@}

  private:

    // Point handle.
    HandleType d_handle;

    // Point coordinates.
    CoordinateType d_x_coord;
    CoordinateType d_y_coord;
    CoordinateType d_z_coord;

  public:
    //! Constructor.
    Point(HandleType _handle,
	  CoordinateType _x, 
	  CoordinateType _y, 
	  CoordinateType _z)
	: d_handle(_handle)
    {   
	d_x_coord = _x;
	d_y_coord = _y;
	d_z_coord = _z;
    }

    //! Copy constructor.
    template<class HandleType_2, class CoordinateType_2>
    Point(const Point<HandleType_2,CoordinateType_2>& point)
        : d_handle( point.handle )
	, d_x_coord( point.d_x_coord )
	, d_y_coord( point.d_y_coord )
	, d_z_coord( point.d_z_coord )
    { /* ... */ }

    //@{
    //! Get the handle.
    HandleType  handle() const { return d_handle; }
    HandleType& handle()       { return d_handle; }
    //@}

    //@{
    //! Get the x coordinate.
    CoordinateType  x() const { return d_x_coord; }
    CoordinateType& x()       { return d_x_coord; }
    //@}

    //@{
    //! Get the y coordinate.
    CoordinateType  y() const { return d_y_coord; }
    CoordinateType& y()       { return d_y_coord; }
    //@}
    
    //@{
    //! Get the z coordinate.
    CoordinateType  z() const { return d_z_coord; }
    CoordinateType& z()       { return d_z_coord; }
    //@}
};

} // end namespace Coupler

#endif // mesh_Point_hpp

//---------------------------------------------------------------------------//
//              end of mesh/Mesh_Point.hpp
//---------------------------------------------------------------------------//
