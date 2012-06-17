//---------------------------------------------------------------------------//
/*!
 * \file DTK_BoundingBox.hpp
 * \author Stuart R. Slattery
 * \brief Bounding box declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOUNDINGBOX_HPP
#define DTK_BOUNDINGBOX_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerializationTraits.hpp>

namespace DataTransferKit
{

class BoundingBox
{

  public:

    // Default constructor.
    BoundingBox();

    // Constructor
    BoundingBox( const double x_min, const double y_min, const double z_min,
		 const double x_max, const double y_max, const double z_max );

    // Destructor.
    ~BoundingBox();

    // Determine if a point is in the box.
    bool pointInBox( double coords[3] ) const;

    // Get the boundaries of the box.
    Teuchos::Tuple<double,6> getBounds() const
    { return Teuchos::tuple( x_min, y_min, z_min, x_max, y_max, z_max ); }

    // Addition assignment operator overload.
    inline BoundingBox& operator +=( const BoundingBox& box );
    
  private:

    // X min.
    double d_x_min;

    // Y min.
    double d_y_min;

    // Z min.
    double d_z_min;

    // X max.
    double d_x_max;

    // Y max.
    double d_y_max;

    // Z max.
    double d_z_max;
};

//---------------------------------------------------------------------------//
// Inline functions.
//---------------------------------------------------------------------------//
/*!
 * \Brief Addition assignment operator overload.
 */
BoundingBox& BoundingBox::operator +=( const BoundingBox& box )
{
    // When we do a global reduction with a bounding box array, we only want
    // one bounding box in each array segment. They will be unique, but this
    // way we will be sure this is true. If the default constructor was called
    // on this box, then set this box's bounds to the other box. Otherwise
    // we'll return this box. 
    if ( d_x_min == 0.0 && d_y_min == 0.0 && d_z_min == 0.0 &&
	 d_x_max == 0.0 && d_y_max == 0.0 && d_z_max == 0.0 )
    {
	d_x_min = box.d_x_min;
	d_y_min = box.d_y_min;
	d_z_min = box.d_z_min;
	d_x_max = box.d_x_max;
	d_y_max = box.d_y_max;
	d_z_max = box.d_z_max;
    }

    return *this;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Serialization traits.
//---------------------------------------------------------------------------//
namespace Teuchos
{
template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::BoundingBox>
    : public DirectSerializationTraits<Ordinal,DataTransferKit::BoundingBox>
{ /* ... */ };
} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end DTK_BOUNDINGBOX_HPP

//---------------------------------------------------------------------------//
// end DTK_BoundingBox.hpp
//---------------------------------------------------------------------------//

