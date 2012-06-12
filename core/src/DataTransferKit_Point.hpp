//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DataTransferKit_Point.hpp
 * \author Stuart R. Slattery
 * \date   Wed May 25 12:23:57 2011
 * \brief  Cartesian Point class definition
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POINT_HPP
#define DTK_POINT_HPP

#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//===========================================================================//
/*!
 * \class Point
 * \brief General point container.
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

    // Index.
    HandleType d_handle;

    // Coordinates.
    CoordinateType d_coords[DIM];

  public:
    
    //! Default constructor.
    Point()
    { /* ... */ }

    //! Copy constructor.
    Point( const Point<DIM,HandleType,CoordinateType>& pt );

    //! Copy constructor.
    Point<DIM,HandleType,CoordinateType>&
    operator=( const Point<DIM,HandleType,CoordinateType>& pt );

    //! Get the handle.
    HandleType getHandle() const 
    { return d_handle; }

    //! Set the handle.
    void setHandle( const HandleType& handle )
    { d_handle = handle; }

    //! Get the coordinates.
    Teuchos::Tuple<CoordinateType,DIM> getCoords() const
    {
	Teuchos::Tuple<CoordinateType,DIM> coords;
	for ( int n = 0; n < DIM; ++n )
	{
	    coords[n] = d_coords[n];
	}
	return coords;
    }

    //! Set the coordinates.
    void setCoords( const CoordinateType coords[DIM] )
    {
	for ( int n = 0; n < DIM; ++n )
	{
	    d_coords[n] = coords[n];
	}
    }
};

//! Create a 1-D point.
template<typename HandleType, typename CoordinateType> inline
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

} // end namespace DataTransferKit

// Create a 1-D point.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Point<1,HandleType,CoordinateType> 
DataTransferKit::point( const HandleType& handle,
			const CoordinateType& x0 )
{
    Point<1,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[1] = { x0 };
    pt.setCoords( coords );
    return pt;
}

// Create a 2-D point.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Point<2,HandleType,CoordinateType> 
DataTransferKit::point( const HandleType& handle,
			const CoordinateType& x0,
			const CoordinateType& x1 )
{
    Point<2,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[2] = { x0, x1 };
    pt.setCoords( coords );
    return pt;
}

// Create a 3-D point.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Point<3,HandleType,CoordinateType> 
DataTransferKit::point( const HandleType& handle,
			const CoordinateType& x0,
			const CoordinateType& x1,
			const CoordinateType& x2 )
{
    Point<3,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[3] = { x0, x1, x2 };
    pt.setCoords( coords );
    return pt;
}

// Create a 4-D point.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Point<4,HandleType,CoordinateType> 
DataTransferKit::point( const HandleType& handle,
			const CoordinateType& x0,
			const CoordinateType& x1, 
			const CoordinateType& x2,
			const CoordinateType& x3 )
{
    Point<4,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[4] = { x0, x1, x2, x3 };
    pt.setCoords( coords );
    return pt;
}

#endif // DTK_POINT_HPP

//---------------------------------------------------------------------------//
//              end of DataTransferKit_Point.hpp
//---------------------------------------------------------------------------//
