//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Coupler_Point.hpp
 * \author Stuart R. Slattery
 * \date   Wed May 25 12:23:57 2011
 * \brief  Cartesian Point class definition
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_POINT_HPP
#define COUPLER_POINT_HPP

#include <Teuchos_Tuple.hpp>

namespace Coupler
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

template<int DIM, typename HandleType=int, typename CoordinateType=double>
class Point 
{
  private:

    // Point handle.
    HandleType d_handle;

    // Point coordinates.
    CoordinateType d_coords[DIM];

  public:
    
    //! Default constructor.
    inline Point();

    //! Copy constructor.
    Point( const Point<DIM,HandleType,CoordinateType>& pt );

    //! Copy constructor.
    Point<DIM,HandleType,CoordinateType>&
    operator=( const Point<DIM,HandleType,CoordinateType>& pt );

    //! Get the handle.
    HandleType getHandle() const 
    { return d_handle; }

    //! Get the coordinates.
    Teuchos::Tuple<CoordinateType,DIM> getCoords() const;
};

//! Create a 1-D point.
template<typename CoordinateType, typename HandleType> inline
Point<1,HandleType,CoordinateType> point( const HandleType& handle,
					  const CoordinateType& x0 );

//! Create a 2-D point.
template<typename HandleType, typename CoordinateType> inline
Point<2,HandleType,CoordinateType> point( const HandleType& handle,
					  const CoordinateType& x0,
					  const CoordinateType& x1 );

//! Create a 3-D point.
template<typename HandleType, typename CoordinateType> inline
Point<3,HandleType,CoordinateType> point( const HandleType& handle,
					  const CoordinateType& x0,
    					  const CoordinateType& x1,
					  const CoordinateType& x2 );

//! Create a 4-D point.
template<typename HandleType, typename CoordinateType> inline
Point<4,HandleType,CoordinateType> point( const HandleType& handle,
					  const CoordinateType& x0,
					  const CoordinateType& x1, 
					  const CoordinateType& x2,
					  const CoordinateType& x3 );

// Copy constructor.
template<int DIM, typename HandleType, typename CoordinateType>
Point<DIM,HandleType,CoordinateType>::Point( 
    const Point<DIM,HandleType,CoordinateType>& pt )
{
    d_handle = pt.d_handle;
    for ( int n = 0; n < DIM; ++n )
    {
	d_coords[n] = pt.d_coords[n];
    }
}

// Copy constructor.
template<int DIM, typename HandleType, typename CoordinateType>
Point<DIM,HandleType,CoordinateType>&
Point<DIM,HandleType,CoordinateType>::operator=( 
    const Point<DIM,HandleType,CoordinateType>& pt )
{
    d_handle = pt.d_handle;
    for ( int n = 0; n < DIM; ++n )
    {
	d_coords[n] = pt.d_coords[n];
    }
    return *this;
}

} // end namespace Coupler

// Create a 1-D point.
template<typename HandleType, typename CoordinateType> inline
Coupler::Point<1,HandleType,CoordinateType> 
Coupler::point( const HandleType& handle,
		const CoordinateType& x0 )
{
    Point<1,HandleType,CoordinateType> pt;
    pt.d_handle = handle;
    pt.d_coords[0] = x0;
}

// Create a 2-D Coupler::point.
template<typename HandleType, typename CoordinateType> inline
Coupler::Point<2,HandleType,CoordinateType> 
Coupler::point( const HandleType& handle,
		const CoordinateType& x0,
		const CoordinateType& x1 )
{
    Point<2,HandleType,CoordinateType> pt;
    pt.d_handle = handle;
    pt.d_coords[0] = x0;
    pt.d_coords[1] = x1;
}

// Create a 3-D Coupler::point.
template<typename HandleType, typename CoordinateType> inline
Coupler::Point<3,HandleType,CoordinateType> 
Coupler::point( const HandleType& handle,
		const CoordinateType& x0,
		const CoordinateType& x1,
		const CoordinateType& x2 )
{
    Point<3,HandleType,CoordinateType> pt;
    pt.d_handle = handle;
    pt.d_coords[0] = x0;
    pt.d_coords[1] = x1;
    pt.d_coords[2] = x2;
}

// Create a 4-D Coupler::point.
template<typename HandleType, typename CoordinateType> inline
Coupler::Point<4,HandleType,CoordinateType> 
Coupler::point( const HandleType& handle,
		const CoordinateType& x0,
		const CoordinateType& x1, 
		const CoordinateType& x2,
		const CoordinateType& x3 )
{
    Point<4,HandleType,CoordinateType> pt;
    pt.d_handle = handle;
    pt.d_coords[0] = x0;
    pt.d_coords[1] = x1;
    pt.d_coords[2] = x2;
    pt.d_coords[3] = x3;
}

#endif // COUPLER_POINT_HPP

//---------------------------------------------------------------------------//
//              end of Coupler_Point.hpp
//---------------------------------------------------------------------------//
