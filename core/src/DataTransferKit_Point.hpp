//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DataTransferKit_Point.hpp
 * \author Stuart R. Slattery
 * \date   Wed May 25 12:23:57 2011
 * \brief  Cartesian Point class definition
 */
//---------------------------------------------------------------------------//

#ifndef DATATRANSFERKIT_POINT_HPP
#define DATATRANSFERKIT_POINT_HPP

namespace DataTransferKit
{
//===========================================================================//
/*!
 * \class Point
 * \brief Cartesian Point class container. Holds (x,y,z)
 */
/*! 
 * \example test/tstPoint.cc
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

} // end namespace DataTransferKit

#endif // DATATRANSFERKIT_POINT_HPP

//---------------------------------------------------------------------------//
//              end of DataTransferKit_Point.hpp
//---------------------------------------------------------------------------//
