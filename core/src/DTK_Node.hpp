//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DTK_Node.hpp
 * \author Stuart R. Slattery
 * \date   Wed May 25 12:23:57 2011
 * \brief  Cartesian Node class definition
 */
//---------------------------------------------------------------------------//

#ifndef DTK_NODE_HPP
#define DTK_NODE_HPP

#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//===========================================================================//
/*!
 * \class Node
 * \brief General node container.
 */
/*! 
 * \example test/tstNode.cc
 *
 * Test of Node.
 */
//===========================================================================//

template<int DIM, typename HandleType=int, typename CoordinateType=double>
class Node 
{
  private:

    // Index.
    HandleType d_handle;

    // Coordinates.
    CoordinateType d_coords[DIM];

  public:
    
    //! Default constructor.
    Node()
    { /* ... */ }

    //! Copy constructor.
    Node( const Node<DIM,HandleType,CoordinateType>& pt );

    //! Copy constructor.
    Node<DIM,HandleType,CoordinateType>&
    operator=( const Node<DIM,HandleType,CoordinateType>& pt );

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

//! Create a 1-D node.
template<typename HandleType, typename CoordinateType> inline
Node<1,HandleType,CoordinateType> node( const HandleType& handle,
					const CoordinateType& x0 );

//! Create a 2-D node.
template<typename HandleType, typename CoordinateType> inline
Node<2,HandleType,CoordinateType> node( const HandleType& handle,
					const CoordinateType& x0,
					const CoordinateType& x1 );

//! Create a 3-D node.
template<typename HandleType, typename CoordinateType> inline
√Node<3,HandleType,CoordinateType> node( const HandleType& handle,
					 const CoordinateType& x0,
					 const CoordinateType& x1,
					 const CoordinateType& x2 );

//! Create a 4-D node.
template<typename HandleType, typename CoordinateType> inline
Node<4,HandleType,CoordinateType> node( const HandleType& handle,
					const CoordinateType& x0,
					const CoordinateType& x1, 
					const CoordinateType& x2,
					const CoordinateType& x3 );

// Copy constructor.
template<int DIM, typename HandleType, typename CoordinateType>
Node<DIM,HandleType,CoordinateType>::Node( 
    const Node<DIM,HandleType,CoordinateType>& pt )
{
    d_handle = pt.d_handle;
    for ( int n = 0; n < DIM; ++n )
    {
	d_coords[n] = pt.d_coords[n];
    }
}

// Copy constructor.
template<int DIM, typename HandleType, typename CoordinateType>
Node<DIM,HandleType,CoordinateType>&
Node<DIM,HandleType,CoordinateType>::operator=( 
    const Node<DIM,HandleType,CoordinateType>& pt )
{
    d_handle = pt.d_handle;
    for ( int n = 0; n < DIM; ++n )
    {
	d_coords[n] = pt.d_coords[n];
    }
    return *this;
}

} // end namespace DataTransferKit

// Create a 1-D node.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Node<1,HandleType,CoordinateType> 
DataTransferKit::node( const HandleType& handle,
		       const CoordinateType& x0 )
{
    Node<1,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[1] = { x0 };
    pt.setCoords( coords );
    return pt;
}

// Create a 2-D node.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Node<2,HandleType,CoordinateType> 
DataTransferKit::node( const HandleType& handle,
		       const CoordinateType& x0,
		       const CoordinateType& x1 )
{
    Node<2,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[2] = { x0, x1 };
    pt.setCoords( coords );
    return pt;
}

// Create a 3-D node.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Node<3,HandleType,CoordinateType> 
DataTransferKit::node( const HandleType& handle,
		       const CoordinateType& x0,
		       const CoordinateType& x1,
		       const CoordinateType& x2 )
{
    Node<3,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[3] = { x0, x1, x2 };
    pt.setCoords( coords );
    return pt;
}

// Create a 4-D node.
template<typename HandleType, typename CoordinateType> inline
DataTransferKit::Node<4,HandleType,CoordinateType> 
DataTransferKit::node( const HandleType& handle,
		       const CoordinateType& x0,
		       const CoordinateType& x1, 
		       const CoordinateType& x2,
		       const CoordinateType& x3 )
{
    Node<4,HandleType,CoordinateType> pt;
    pt.setHandle( handle );
    CoordinateType coords[4] = { x0, x1, x2, x3 };
    pt.setCoords( coords );
    return pt;
}

#endif // DTK_NODE_HPP

//---------------------------------------------------------------------------//
//              end of DTK_Node.hpp
//---------------------------------------------------------------------------//
